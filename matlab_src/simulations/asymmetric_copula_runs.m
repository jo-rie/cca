function asymmetric_copula_runs(exp_path, ns, nb_values_to_save, a_values, b_values)
    

    copula_cdf = @(a, b, u, v) u .* v + u .* v .* (1-u) .* (1-v) .* ((a-b) .* v .* (1-u) + b);

    % b = .5; 
    % cond_on_a = @(b) b - 3 - (9 + 6 .* b - 3 .* b .^2).^(1/2) / 2;
    % plot(-1:0.01:1, cond_on_a(-1:0.01:1))
    % fprintf('cond: a ≥ %.4f\n', b - 3 - (9 + 6 * b - 3 * b .^2)^(1/2) / 2)

    for n = ns
        sing_vecs_left = zeros(numel(a_values), n, nb_values_to_save);
        sing_vals = zeros(numel(a_values), nb_values_to_save);
        sing_vecs_right = zeros(numel(a_values), n, nb_values_to_save);
        rowProfiles = zeros(numel(a_values), n, 2);
        columnProfiles = zeros(numel(a_values), n, 2);
        for i = 1:2
            a = a_values(i);
            b = b_values(i);
            disc_cop = CopulaCAByCdfFunction(@(u, v) copula_cdf(a, b, u, v), n);
            sing_vecs_left(i, :, :) = getU(disc_cop, nb_values_to_save);
            sing_vals(i, :) = getSingVals(disc_cop, nb_values_to_save);
            sing_vecs_right(i, :, :) = getV(disc_cop, nb_values_to_save);
            rowProfiles(i, :, :) = disc_cop.getRowProfiles(2);
            columnProfiles(i, :, :) = disc_cop.getColumnProfiles(2);
            plot_grid = disc_cop.plotGrid;
        end
    end

    save(fullfile(exp_path, sprintf("results_n_%i.mat", n)), ...
                'sing_vals', "sing_vecs_left", "sing_vecs_right", "a_values", ...
                "b_values", "n", "plot_grid", 'nb_values_to_save', 'rowProfiles', 'columnProfiles'); 
end % asymmetric_copula_runs