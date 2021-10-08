function [v_out, cost, times, w_out] = isospec(M, N, par, w_init)

manifold = euclideanfactory(M.n, 1);
problem = {};
problem.M = manifold;
problem.cost = @(v) obj_std_tanh2(v, M, N, par.k, par.tau);
problem.egrad = @(v) grad_std_tanh2(v, M, N, par.k, par.tau);
% figure, checkgradient(problem)
    
options.tolgradnorm = 1e-10;
options.maxiter = par.iter;
options.maxtime = 5 * 60;
options.verbosity = 1;

w_out = zeros(M.n, par.n_start);

dV = par.tau;

tic
parfor i=1:par.n_start
    
    w_out(:, i) = trustregions(problem, w_init(:, i), options);
    
    v_out(:, i) = 0.5 * dV * (tanh(w_out(:, i)) + 1);
    
    cost(i) = problem.cost(w_out(:, i));
end
times = toc;

end