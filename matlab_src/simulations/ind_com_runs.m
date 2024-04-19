function ind_com_runs(exp_path, n)

    comon_cdf = @(u, v) min(u, v);

    comon_cop = CopulaCAByCdfFunction(comon_cdf, n);
    rowprofiles = comon_cop.getRowProfiles(2);
    k = 1;
    [eta, fro_er] = mar_minimize_eta(comon_cop.getPDFMatrix(), k);
    comon_mar_cop = CopulaMAR(comon_cdf, n, eta);
    rowprofiles_MAR = comon_mar_cop.getRowProfiles(2);

    plotGrid = comon_cop.plotGrid;

    save(fullfile(exp_path, 'com_mar.mat'), 'rowprofiles', 'rowprofiles_MAR', 'plotGrid', 'k', 'n', 'eta');

    ind_cdf = @(u, v) u.* v;
    ind_cop = CopulaCAByCdfFunction(ind_cdf, n);
    rowprofiles = ind_cop.getRowProfiles(2);
    plotGrid = ind_cop.plotGrid;

    save(fullfile(exp_path, 'ind.mat'), 'rowprofiles', 'n', 'plotGrid');

end % ind_com_runs