function [fro] = mar_fro4eta(eta, k, p)
    %% MAR_FRO4ETA Compute the resulting frobenius norm between the pdf-matrix p and the MAR model for parameter eta for truncation of rank k

    n = length(p);
    func_p_eta = @(eta) p + ...
        eta * eye(n) - ...
        (1 + eta) ./ n .* ones(n);
    func_p_eta_inv = @(p_ab, eta) (p_ab - eta .* eye(n) + (1 + eta) ./ n .* ones(n));

    p_eta = func_p_eta(eta);
    [U, S, V] = svd(p_eta);
    p_eta_rec = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
    p_rec = func_p_eta_inv(p_eta_rec, eta);
%     p_rec = a 
    fro = norm(p - p_rec, 'fro');
end