function mask = grow(M, seed, n_edges)

if n_edges<0
    error('Cannot grow a negative number of steps.')
end

onering = calc_onering(M);

%FIXME: probably i don't need dilated as a matrix, just as a vector

dilated = false(M.n,n_edges);
mask = seed;

if n_edges==0
    return;
end

for k=1:n_edges
    for i=1:M.n
        if ~mask(i) && ~any(dilated(i,1:k-1))
            if any(mask(onering{i})) || any(any(dilated(onering{i},1:k-1),2))
                dilated(i,k) = 1;
            end
        end
    end
    mask = mask | dilated(:,k);
end

end
