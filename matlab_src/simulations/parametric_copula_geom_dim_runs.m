function parametric_copula_geom_dim_runs(exp_path, copulas, theta_ranges, ns)

    for i_copulas = 1:length(copulas)
        fprintf('Starting %s\n', copulas(i_copulas));
        exp_path_cop = fullfile(exp_path, copulas(i_copulas));
        [~,~,~] = mkdir(exp_path_cop);
    
        thetas = theta_ranges{i_copulas};
        %%%%%% Geometric Dimension computation
        geom_dim = zeros(length(ns), length(thetas));
%        fprintf("-> Starting geometric dimension computations for %s\n", exp_path)
        textprogressbar('Running...')
    
        for i_n = 1:length(ns)
            for i_theta = 1:length(thetas)
                textprogressbar(((i_n - 1) * length(thetas) + i_theta - 1) / length(ns) / length(thetas) * 100)
                n = ns(i_n);
                theta = thetas(i_theta);
                disc_cop = CopulaCAParametric(copulas(i_copulas), theta, n, 'theta');
                geom_dim(i_n, i_theta) = disc_cop.getGeometricDimension();
            end
        end
    
        textprogressbar(' Done.')
        save(fullfile(exp_path_cop, 'geom_dim_n_theta.mat'), 'ns', 'thetas', 'geom_dim');


    end % i_copulas
end % parametric_copula_geom_dim_runs