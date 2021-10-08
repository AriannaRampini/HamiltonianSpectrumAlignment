function x = obj_std_tanh2(w, M, N, k, prof_buca)

e = 1;
try
    evals = eigs(M.W + M.A*spdiag(0.5 * prof_buca * (tanh(w) + 1)), M.A, k, 'SM');
catch
    fprintf("cond(M.W) = %f\n cond(pot) = %f\n cond(somma) = %f\n", cond(M.W, 1), cond(M.A*spdiag(0.5 * prof_buca * (tanh(w) + 1)), 1), ...
        cond(M.W + M.A*spdiag(0.5 * prof_buca * (tanh(w) + 1)), 1));
    evals = eigs(M.W + M.A*spdiag(0.5 * prof_buca * (tanh(w) + 1)) + spdiag(e * ones(M.n,1)), M.A, k, 'SM');
end
a = (N.evals(1:k).^2);
x = sum( ((evals - N.evals(1:k)).^2)./a); 

end
