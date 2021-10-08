function mask = baseline_kmeans(M, N)

G = calc_shot(M.VERT', M.TRIV', 1:M.n, 9, 10, 3)';
G(:, ~any(G,1)) = [];  % Delete all zero descriptors

F = calc_shot(N.VERT', N.TRIV', 1:N.n, 9, 10, 3)';
F(:, ~any(F,1)) = [];  % Delete all zero descriptors

tmp = kmeans([G; F], min(4000, N.n));
mask = ismember(tmp(1:M.n), unique(tmp((M.n+1):(M.n+N.n))));

end