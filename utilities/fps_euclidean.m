function [ S ] = fps_euclidean(V, n, seed)
%EUCLIDEANFPS Samples K vertices from V by using farthest point sampling.
% The farthest point sampling starts with vertex seed and uses the euclidean
% metric of the 3-dimensional embedding space.
% -  V is a n-by-3 matrix storing the positions of n vertices
% -  K is the number of samples
% -  seed is the index of the first vertex in the sample set
% Returns
% -  S is a K-dimensional vector that includes the indices of the K sample
%    vertices.

S = zeros(n,1);
S(1) = seed;
d = pdist2(V,V(seed,:));

for i=2:n
    [~,m] = max(d);
    S(i) = m(1);
    d = min(pdist2(V,V(S(i),:)) , d);
end

end
