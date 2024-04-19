function frobenius_different_families(exp_path, cops, cop_short, fro_n)

results_fro = zeros(length(cops));
results_d1 = zeros(length(cops));
results_d2 = zeros(length(cops));

for i = 1:length(cops)
    for j = 1:length(cops)
        cop1 = cops{i};
        cop2 = cops{j};
        val_fro = dfrobenius(cop1, cop2, fro_n);
        results_fro(i, j) = val_fro;
        results_fro(j, i) = val_fro;
        val_d1 = d1(cop1, cop2, fro_n);
        results_d1(i, j) = val_d1;
        results_d1(j, i) = val_d1;
        val_d2 = d2(cop1, cop2, fro_n);
        results_d2(i, j) = val_d2;
        results_d2(j, i) = val_d2;
    end
end
save(fullfile(exp_path, sprintf('frobenius_different_families_n_%i.mat', fro_n)), '-mat', ...
        'fro_n', 'cop_short', 'results_fro', ...
        'results_d2', 'results_d1')
end