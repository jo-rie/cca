function frobenius_gauss_runs(exp_path, rhos, n)
    results_fro = zeros(length(rhos));
    results_d1 = zeros(length(rhos));
    results_d2 = zeros(length(rhos));
    for i = 1:length(rhos)
        for j = 1:length(rhos)
            if i > j
                cop1 = CopulaCAParametric('Gaussian', rhos(i), n, 'theta');
                cop2 = CopulaCAParametric('Gaussian', rhos(j), n, 'theta');
                val_fro = dfrobenius(cop1, cop2, n);
                results_fro(i, j) = val_fro;
                results_fro(j, i) = val_fro;
                val_d1 = d1(cop1, cop2, n);
                results_d1(i, j) = val_d1;
                results_d1(j, i) = val_d1;
                val_d2 = d2(cop1, cop2, n);
                results_d2(i, j) = val_d2;
                results_d2(j, i) = val_d2;
            end
        end
    end
    save(fullfile(exp_path, sprintf('distance_n_%i_rhos.mat', n)), '-mat', ...
        'n', 'rhos', 'results_fro', 'results_d1', 'results_d2')

end %frobenius_gauss_runs
