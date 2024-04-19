function hellinger_run(cop_fam, cop_theta, cop_n, exp_path)
%HELLINGER_RUN Hellinger comparison for specified copula

cop_obj = CopulaCAParametric(cop_fam, cop_theta, cop_n, 'theta');
pdf_matrix = cop_obj.getPDFMatrix();
plot_grid = cop_obj.plotGrid;

% Compute Parameters for Frobenius

sing_vecs_fro = cop_obj.getU(5);
sing_vals_fro = cop_obj.getSingVals(5);

% Compute Parameters for Hellinger

[U, S, V] = svds(pdf_matrix .^ 0.5, 5);
for i = 1:5
    if U(1, i) < 0
        U(:, i) = - U(:, i);
    end
end
sing_vecs_hellinger = U;
sing_vals_hellinger = diag(S);

checkerboard_fro = cop_obj.getApproxPDFMatrix(5);
checkerboard_hellinger = (U * S * V') .^ 2;

save(fullfile(exp_path, sprintf('%s_theta_%.2f_n_%i.mat', ...
    cop_fam, cop_theta, cop_n)), ...
    "checkerboard_hellinger", "checkerboard_fro", "sing_vals_hellinger", ...
    "sing_vecs_hellinger", "sing_vals_fro", "sing_vecs_fro", 'plot_grid')

end

