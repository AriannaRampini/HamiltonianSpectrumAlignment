function dv = grad_std_tanh2(w, M, N, k, prof_buca)

[evecs, evals] = eigs(M.W + M.A*spdiag(0.5 * prof_buca * (tanh(w) + 1)), M.A, k, 'SM');
Q = evals - spdiag(N.evals(1:k));
a = (N.evals(1:k).^2)'; 
dv = prof_buca * (1 - (tanh(w)).^2) .* sum(((M.A*(evecs.^2))*Q)./a, 2);

end
