function frobenius_over_k(exp_path, tau1, tau2, n)
    val_d1 = zeros(length(tau1), n);
    val_d2 = zeros(length(tau1), n);
    for i = 1:length(tau1)
        cop1 = CopulaCAParametric('Gaussian', tau1(i), n);
        cop2 = CopulaCAParametric('Gaussian', tau2(i), n);
        for k = 1:n
            val_d1(i, k) = d1(cop1, cop2, k);
            val_d2(i, k) = d2(cop1, cop2, k);
        end
    end
    
    save(fullfile(exp_path, 'Gaussian_over_k.mat'), '-mat', 'tau1', 'tau2', 'n', 'val_d1', 'val_d2');
end