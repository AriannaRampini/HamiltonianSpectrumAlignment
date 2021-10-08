function init = multistart(M, n, sigma, seed)

init = ones(M.n, n);
centers = fps_euclidean(M.VERT, n, seed);

for i=1:n
    init(:, i) = init(:, i) - exp(-sum((M.VERT - M.VERT(centers(i), :)).^2, 2)/2/sigma^2);
end

end