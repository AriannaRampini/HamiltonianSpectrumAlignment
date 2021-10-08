function nn = nearest_neig(M, centers, order)
    % trova sulla mesh M i vertici più vicini ai centri centers nello
    % spazio euclideo. se order, mi dice anche a quale centro corrisponde 
    % il vertice 

    nn = zeros(M.n, 1);
    for i = 1:size(centers, 1)
        [~, idx_min] = min(sum((M.VERT - centers(i, :)).^2, 2));
        nn(idx_min) = nn(idx_min) + 1;
    end
    if not(order)
        nn = logical(nn);
    end
    
end