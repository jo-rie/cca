function parametric_copula_mar_runs(exp_path, copulas, theta_ranges, ns, ...
    approx_orders)

    for i_copulas = 1:length(copulas)
        fprintf('Starting %s\n', copulas(i_copulas));
        exp_path_cop = fullfile(exp_path, lower(copulas(i_copulas)));
        [~,~,~] = mkdir(exp_path_cop);
    
        thetas = theta_ranges{i_copulas};
        %%%%%% Computation of eta values
        fprintf("-> Starting MAR computations\n")    
        for n = ns
            textprogressbar(sprintf('-> Running n = %i: ', n))
            etas = zeros(length(thetas), length(approx_orders));
            frobenius_error_mar = zeros(length(thetas), length(approx_orders));
            frobenius_error_raw = zeros(length(thetas), length(approx_orders));
            rowprofiles = zeros(length(thetas), n, 2);
            rowprofiles_MAR = zeros(length(thetas), n, 2);
            for i_theta = 1:length(thetas)
                for appr_ord_i = 1:length(approx_orders) 
                    theta = thetas(i_theta);
                    approx_order = approx_orders(appr_ord_i);
                    textprogressbar(((i_theta-1) * length(approx_orders) + appr_ord_i - 1) / length(thetas) / length(approx_orders) * 100) 
                    disc_cop = CopulaCAParametric(copulas(i_copulas), theta, n, 'theta');
                    [eta, fro_er] = mar_minimize_eta(disc_cop.getPDFMatrix(), approx_order); 
                    fro_er_raw = disc_cop.getPDFMatrixApproxError(approx_order, false).('fro');
                    mar_cop = CopulaMAR(@(u, v) mycopulacdf(copulas(i_copulas), [u, v], theta), n, eta);
            %         fprintf(['Theta: %.2f, a = %.4f, b = %.4f, Fro-Error = %.4f, ' ...
            %             'Fro-Error raw = %.4f\n'], theta, alpha, beta, fro_er, ...
            %             fro_er_raw)
                    etas(i_theta, appr_ord_i) = eta;
                    frobenius_error_mar(i_theta, appr_ord_i) = fro_er;
                    frobenius_error_raw(i_theta, appr_ord_i) = fro_er_raw;
                    rowprofiles(i_theta, :, :) = disc_cop.getRowProfiles(2);
                    rowprofiles_MAR(i_theta, :, :) = mar_cop.getRowProfiles(2);
                end
            end
            save(fullfile(exp_path_cop, sprintf('mar_n_%i.mat', n)), 'etas', ...
                'frobenius_error_raw', 'frobenius_error_mar', 'thetas', ...
                'n', 'approx_orders', 'rowprofiles', 'rowprofiles_MAR')
            textprogressbar(' Done.')
        end
    end
end
