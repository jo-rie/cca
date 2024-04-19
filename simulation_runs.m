export_root = fullfile(config_local.base_path, ...
    config.results_dir);
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'stats:copulastat:BadGumbelParameter');
warning('on', 'stats:internal:getParamVal:BadValue');

%% Analysis of Independence and Comonotonicity Copula

ind_com_exp_path = fullfile(export_root, 'ind_com');
[~,~,~] = mkdir(ind_com_exp_path);

ind_com_n = 50;

ind_com_runs(ind_com_exp_path, ind_com_n);


clear -regexp ^ind_com;

%% Asymmetric Copula Runs

asym_ns = [50, 100];
asym_nb_values_to_save = 10;
asym_a_values = [-1.5, .5];
asym_b_values = [.5, -.5];
asym_exp_path = fullfile(export_root, 'asymmetric_copula');
[~] = mkdir(asym_exp_path);

asymmetric_copula_runs(asym_exp_path, asym_ns, asym_nb_values_to_save, ...
    asym_a_values, asym_b_values)
clear -regexp ^asym;

%% Parametric Copula Runs

param_copulas = ["ca", "amh", "gumbel", "gaussian", "fgm"];
param_theta_ranges_full = {...
        0.05:0.05:0.95, ...
        -1:0.1:1, ...
        1:0.5:10, ...
        -0.75:0.25:0.75, ...
        -1:0.05:1
    };
param_theta_ranges_geom_dim = {...
    [0.25, 0.5, 0.75],...  % ca
    [0.3, 0.5, 0.7, 0.9],...  % amh
    [2.5, 5, 7.5, 10],...  % gumbel
    [0.25, 0.5, 0.75],...  % gaussian
    [.8]...  % fgm
};
param_ns = [50, 100];
param_ns_geom_dim = round(logspace(1, log10(10000), 20));
param_nb_values_to_save = 10;
param_approx_orders = 0:5;
param_exp_path = export_root;
parametric_copula_raw_runs(param_exp_path, param_copulas, ...
    param_theta_ranges_full, param_ns, param_nb_values_to_save)

fprintf('Starting MAR runs at %s \n', datetime('now'))
parametric_copula_mar_runs(param_exp_path, param_copulas, param_theta_ranges_full, param_ns, param_approx_orders)

fprintf('Starting geometric dimension runs at %s \n', datetime('now'))
parametric_copula_geom_dim_runs(param_exp_path, param_copulas, param_theta_ranges_geom_dim, param_ns_geom_dim)

% Separate runs for FGM copula with high n
fprintf('Starting separate FGM  at %s \n', datetime('now'))
parametric_copula_raw_runs( ...
    param_exp_path, ...
    ["fgm"], ...
    param_theta_ranges_geom_dim(5), ...
    [1000, 10000], ...
    param_nb_values_to_save...
    )

% Separate runs for Gaussian copula with high n
fprintf('Starting Gaussian copula runs at %s \n', datetime('now'))
gaussian_sing_vecs_for_high_n(...
    [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], ...
    cell2mat(param_theta_ranges_geom_dim(4)),...
    8, ...
    fullfile(export_root, 'gaussian')...
    );
clear -regexp ^param;

%% Correspondence Analysis Plots

cora_n = 100;
cora_param_copulas = ["AMH", "Gumbel", "Gaussian", "Clayton"];
cora_param_theta_ranges = {...
        [-.15, -.05, .05, .25], ...
        .1:.2:.9, ...
        .1:.2:.9, ...
        .1:.2:.9 ...
    };
cora_exp_path = fullfile(export_root, "correspondence_analysis_plots");


clear -regexp ^cora

%% Counterexample Geometric Dimension

cg_center_points = ([...
    0, 0;... %1
    1, 1;...
    4, 1;...
    5, 0;...
    3, 2;... %5
    6, 2;...
    2, 3;...
    7, 3;...
    1, 4;... %9
    4, 4;...
    0, 5;...
    5, 5;...
    2, 6;... %13
    7, 6;...
    3, 7;...
    6, 7
].* 2 + 1)./16;

% Create list of polygons
cg_polyshapes = repmat(polyshape, length(cg_center_points), 1);
cg_area_sum = 0;
for k = 1:length(cg_center_points)
    cg_polyshapes(k) = nsidedpoly(4, 'Center', cg_center_points(k, :), 'SideLength', 1/8);
    % fprintf('Area %i: %.4f\n', k, area(cg_polyshapes(k)))
    cg_area_sum = cg_area_sum + area(cg_polyshapes(k));
end

cg_ns = 1:1:8;
cg_exp_path = fullfile(export_root, 'counterexample_geometric_dimension');
[~, ~, ~] = mkdir(cg_exp_path);

counterexample_polygon_runs(cg_exp_path, cg_center_points, cg_polyshapes, cg_ns)

clear -regexp ^cg;

%% Plots Frobenius Distance

fro_rhos_gauss = (0:0.05:.95) + 0.025;
fro_n = 100;
fro_exp_path = fullfile(export_root, 'frobenius_distance');
[~, ~, ~] = mkdir(fro_exp_path);

frobenius_gauss_runs(fro_exp_path, fro_rhos_gauss, fro_n)

fro_tau1 = [.1, .1, .3, .3]; % Build up pairs of tau
fro_tau2 = [.9, .7, .7, .5];

frobenius_over_k(fro_exp_path, fro_tau1, fro_tau2, fro_n);

fro_cop_short = ["G, 0.3", "AMH, 0.3", "Cl, 0.3", "Ga, 0.3", "G, 0.9", "Cl, 0.9", "Ga, 0.9"];
fro_tau1 = .3;
fro_tau2 = .9;

fro_cops = {CopulaCAParametric('AMH', fro_tau1, fro_n), ...
    CopulaCAParametric("Gumbel", fro_tau1, fro_n), ...
    CopulaCAParametric('Clayton', fro_tau1, fro_n), ...
    CopulaCAParametric('Gaussian', fro_tau1, fro_n), ...
    CopulaCAParametric("Gumbel", fro_tau2, fro_n), ...
    CopulaCAParametric('Clayton', fro_tau2, fro_n), ...
    CopulaCAParametric('Gaussian', fro_tau2, fro_n)};

frobenius_different_families(fro_exp_path, fro_cops, fro_cop_short, fro_n)

frobenius_increasing_n(fro_exp_path, round(logspace(1, 4, 50)))

clear -regexp ^fro

%% Plot Invalid Copula

ic_family = ["Gumbel", "Gumbel"];
ic_theta = [2.5, 7.5];
ic_n = [50, 50];
ic_k = [10, 10];

ic_exp_path = fullfile(export_root, 'invalid_values');
[~, ~, ~] = mkdir(ic_exp_path);

for i = 1:2
    invalid_copula_vals(...
        ic_exp_path, ic_family(i), ic_theta(i), ic_n(i), ic_k(i));
end

clear -regexp ^ic
%% Hellinger analysis

ha_exp_path = fullfile(export_root, 'hellinger');
[~,~,~] = mkdir(ha_exp_path);
ha_cop = "Gaussian";
ha_rho = 0.75;
ha_n = 100;
hellinger_run(ha_cop, ha_rho, ha_n, ha_exp_path);

clear -regexp ^ha





