function parametric_copula_raw_runs(exp_path, copulas, theta_ranges, ns, ...
    nb_values_to_save)

    for i_copulas = 1:length(copulas)
        fprintf('Starting %s\n', copulas(i_copulas));
        exp_path_cop = fullfile(exp_path, lower(copulas(i_copulas)));
        [~] = mkdir(exp_path_cop);
    
        thetas = theta_ranges{i_copulas};
    
        %%%%%% Compute singular vectors and singular values
        for n = ns
            sing_vecs = zeros(length(thetas), n, nb_values_to_save);
            sing_vals = zeros(length(thetas), nb_values_to_save);
            rowprofiles = zeros(length(thetas), n, 2);
            textprogressbar(sprintf('-> Running n = %i: ', n))
            for i = 1:length(thetas)
                theta = thetas(i);
                textprogressbar(i / length(thetas) * 100)
                disc_cop = CopulaCAParametric(copulas(i_copulas), theta, n, 'theta');
                sing_vecs(i, :, :) = getU(disc_cop, nb_values_to_save);
                sing_vals(i, :) = getSingVals(disc_cop, nb_values_to_save);
                rowprofiles(i, :, :) = disc_cop.getRowProfiles(2);
                plot_grid = disc_cop.plotGrid;
                       
            end
            save(fullfile(exp_path_cop, sprintf("raw_n_%i.mat", n)), ...
                'sing_vals', "sing_vecs", "thetas", "n", "plot_grid", ...
                'nb_values_to_save', "rowprofiles");
            textprogressbar('Done')
        end
    end
end