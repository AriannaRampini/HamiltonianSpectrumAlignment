%%
clc
clear all
close all

addpath(genpath('./utilities'))
addpath(genpath('./manopt'))

%% Set parameters

step_coeff = 50;         % coefficient for computing step potential height

par1.k = 20;             % number of eigenvalues used
par1.n_start = 41;       % number of initializations used
par1.iter = 100;         % number of iteration

k_desc = 5;              % number of eigenfunctions used to choose the final solution
frac_vertices = 0.5;     % fraction of vertices kept in the template

%% Load data

M = load_off('./data/cat.off');                % template
N = load_off('./data/cuts_cat.off');           % partial shape

plot_input(M, N)

%% Compute Laplacian

[M.W,~,M.A] = calc_LB_FEM_bc(M, 'dirichlet');
                
[N.W,~,N.A] = calc_LB_FEM_bc(N, 'dirichlet');
[N.evecs, N.evals] = eigs(N.W, N.A, par1.k, 'SM');
N.evals = diag(N.evals);
[N.evals, idx] = sort(N.evals);
N.evecs = N.evecs(:,idx);

%% Compute shot descriptors

M.area = sum(calc_tri_areas(M));
SHOT_M = calc_shot(M.VERT', M.TRIV', 1:M.n, 9, 10, 3)'; 
SHOT_M(:, ~any(SHOT_M,1)) = [];

N.area = sum(calc_tri_areas(N));
SHOT_N = calc_shot(N.VERT', N.TRIV', 1:N.n, 9, 10, 3)';
SHOT_N(:, ~any(SHOT_N,1)) = [];

%% Remesh template

n_remesh = floor(M.n * frac_vertices);
grow_steps = 4;
Mr = remesh(M, struct('vertices', n_remesh));
[Mr.W,~,Mr.A] = calc_LB_FEM_bc(Mr, 'dirichlet');
    
%% Initialization

% Compute the baseline solution and project it on the remeshed template:

R_found_kmeans = baseline_kmeans(M, N);
init_base_rem = nearest_neig(Mr,  M.VERT(logical(R_found_kmeans), :), 0);

% Gaussians centered around farthest samples on M, with 2 difference variance:

init1 = multistart(Mr, (par1.n_start - 1)/4, sqrt(N.area), 1);
init2 = multistart(Mr, (par1.n_start - 1)/4, 2*sqrt(N.area), 1);

w_init = [(2 - 4*init_base_rem), ...
    (init1 - 0.5), (init2 - 0.5), (-init1 + 0.5), (-init2 + 0.5)];


%% Optimization

par1.tau = step_coeff * N.evals(par1.k);

[v_out, cost, times, w_out] = isospec(Mr, N, par1, w_init);
                
%% Choose final solution

% projection of shot descriptors of N on the first k_desc Dirichlet 
% eigenfunctions:

coeffSHOT_N = (N.evecs(:, 1:k_desc).^2)' * N.A * SHOT_N;

% projection of shot descriptors of M on the first k_desc Hamiltonian
% eigenfunctions:

score_SHOT_coeff = ones(1, par1.n_start);
R_found = zeros(M.n, par1.n_start);

parfor i = 1:par1.n_start
    [evecs_out, evals_out] = eigs(Mr.W + Mr.A*spdiag(v_out(:, i)), Mr.A, par1.k, 'SM');
    evals_out = diag(evals_out);
    [~, idx] = sort(evals_out);
    evecs_out = evecs_out(:, idx);
    
    R_found_rem = sum(evecs_out.^2, 2)/par1.k > 0.01*evals_out(1);
    nn = nearest_neig(M,  Mr.VERT(R_found_rem, :), 0);
    R_found(:, i) = grow(M, nn, grow_steps);
    R_found(:, i) = shrink(M, R_found(:, i), grow_steps);
    
    [evecs_out, evals_out] = eigs(M.W + M.A*spdiag(par1.tau * (1 - R_found(:, i))), M.A, k_desc, 'SM');
    evals_out = diag(evals_out);
    [~, idx] = sort(evals_out);
    evecs_out = evecs_out(:,idx);
    coeffSHOT_M = (evecs_out(:, 1:k_desc).^2)' * M.A * SHOT_M;
    score_SHOT_coeff(i) = sum(sqrt( sum((coeffSHOT_M - coeffSHOT_N).^2, 1) ), 2);
end

[~, idx_shot_coeff] = min(score_SHOT_coeff, [], 2);

R_found_final = R_found(:, idx_shot_coeff);

%% Plot result

plot_mask(M, N, R_found_final)
% plot_mask(M, N, R_found_kmeans)
