function frobenius_increasing_n(exp_path,n_array)
% Function to compute for the same pair of copulas with increasing n the
% Frobenius distance and d1 and d2 distances
cop1_cellarray_from_n = {
    @(n) CopulaCAParametric('Gaussian', 0.1, n, 'theta'), ...
    @(n) CopulaCAParametric('Gaussian', 0.4, n, 'theta'), ...
    @(n) CopulaCAParametric('Gumbel', 1, n, 'theta'), ...
    @(n) CopulaCAParametric('Gumbel', 4, n, 'theta')
};

cop2_cellarray_from_n = {
    @(n) CopulaCAParametric('Gaussian', 0.9, n, 'theta'), ...
    @(n) CopulaCAParametric('Gaussian', 0.6, n, 'theta'), ...
    @(n) CopulaCAParametric('Gumbel', 10, n, 'theta'), ...
    @(n) CopulaCAParametric('Gumbel', 6, n, 'theta')
};

cop_str = { ...
    "Gaussian0.1;Gaussian0.9",...
    "Gaussian0.4,Gaussian0.6",...
    "Gumbel1;Gumbel10", ...
    "Gumbel4;Gumbel6"...
    };

fro_er_array = zeros(length(cop1_cellarray_from_n), length(n_array));
d1_array = zeros(length(cop1_cellarray_from_n), length(n_array));
d2_array = zeros(length(cop1_cellarray_from_n), length(n_array));

for i_cops = 1:length(cop1_cellarray_from_n)
    fprintf('Starting %s, n = ', cop_str{i_cops})
    cop1_from_n = cop1_cellarray_from_n{i_cops};
    cop2_from_n = cop2_cellarray_from_n{i_cops};

    for i_n = 1:length(n_array)
        n = n_array(i_n);
        fprintf('%i ...', n)
        cop1 = cop1_from_n(n);
        cop2 = cop2_from_n(n);
        fro_er_array(i_cops, i_n) = dfrobenius(cop1, cop2, n);
        d1_array(i_cops, i_n) = d1(cop1, cop2, n);
        d2_array(i_cops, i_n) = d2(cop1, cop2, n);
    end
    fprintf('\n')
end


save(fullfile(exp_path, 'frobenius_increasing_n.mat'), 'fro_er_array', "d1_array", "d2_array", "n_array", "cop_str")
end

