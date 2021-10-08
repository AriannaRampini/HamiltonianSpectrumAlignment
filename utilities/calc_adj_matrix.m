function A = calc_adj_matrix(M, sym)

% asymmetric
A = sparse(...
    [M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,3)], ...
    [M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,1)], ...
    1, ...
    M.n, M.n, 3 * M.m);

if nargin==2 && sym
    A = A+A';
    A = double(A~=0);
end

end
