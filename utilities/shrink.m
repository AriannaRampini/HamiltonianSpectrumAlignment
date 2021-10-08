function R = shrink(M, R_in, n_step)

    vicini = calc_onering(M);
    R = R_in;
    
    for j = 1:n_step
        flag = zeros(M.n, 1);
        idx_region = find(R);
        
        for i = 1:sum(R)           
            if not(all(R(vicini{idx_region(i)})))
                flag(idx_region(i)) = 1;
            end
        end
        
        R = R - flag;
    end
    
end