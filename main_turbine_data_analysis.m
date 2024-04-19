% Load data

% Setup combinations of time, pressure and size

times = [50, 60, 70, 70, 70, 70, 90, 90, 90, 90];
pressures = [1, 1, 1, 1, 5, 5, 1, 1, 5, 5];
sizes = [640, 640, 230, 640, 230, 640, 230, 640, 230, 640];
base_path = config_local.turbine_base_path;
import_path = fullfile(base_path, 'raw_data');
data_save_path = fullfile(base_path, 'results_copca_20231018');
[~, ~, ~] = mkdir(data_save_path);

grid_size = 101;

%% Iterate over scenarios and compute profiles
for scenario_number = 1:length(times)

% Load data
loop_time = times(scenario_number);
loop_pressure = pressures(scenario_number);
loop_size = sizes(scenario_number);
data = load_file(loop_time, ...
    loop_pressure, loop_size, import_path);

% data_uv = pobs(data); % Is done by EmpiricalCA

loop_export_path = fullfile( ...
    data_save_path, ...
    sprintf('t_%i_p_%i_s_%i_n_%i', times(scenario_number), ...
    pressures(scenario_number), sizes(scenario_number), grid_size));
[~, ~, ~] = mkdir(loop_export_path);


% Iterate over pairs of variables
for i = 1:size(data, 2)
    for j = 1:size(data, 2)
        if i < j
            data_loop = data(:, [i, j]);
            copca_model = CopulaCAEmpirical(data_loop, grid_size);
        
            rowProfiles = copca_model.getRowProfiles(10);
            columnProfiles = copca_model.getColumnProfiles(10);
            num_obs = size(data, 1);
            data_ranked = copca_model.data_ranked;
            checkerboard_pdf = copca_model.getPDFMatrix();

            save_file_name = fullfile(loop_export_path, ...
                sprintf('dim1_%i_dim2_%i.mat', i, j));
            
            save(save_file_name, 'rowProfiles', 'columnProfiles', ...
                "num_obs", 'loop_time', 'loop_pressure', 'loop_size', ...
                "grid_size", 'data_ranked', 'data_loop', 'checkerboard_pdf');
        end %if
    end %j
end %i

fprintf('Finished scenario %i of %i\n', scenario_number, length(times))
end % scenario_number


%% Sample row and column profiles of independence copula

% Check data length for scenarios
% for scenario_number = 1:length(times)
% 
%     % Load data
%     loop_time = times(scenario_number);
%     loop_pressure = pressures(scenario_number);
%     loop_size = sizes(scenario_number);
%     data = load_file(loop_time, ...
%         loop_pressure, loop_size, import_path);
%     fprintf('t_%i_p_%i_s_%i_n_%i: %i\n', times(scenario_number), ...
%     pressures(scenario_number), sizes(scenario_number), grid_size, size(data, 1));
% end % scenario_number

% t_50_p_1_s_640_n_50: 1242
% t_60_p_1_s_640_n_50: 1169
% t_70_p_1_s_230_n_50: 3000
% t_70_p_1_s_640_n_50: 2784
% t_70_p_5_s_230_n_50: 5231
% t_70_p_5_s_640_n_50: 2619
% t_90_p_1_s_230_n_50: 3631
% t_90_p_1_s_640_n_50: 2191
% t_90_p_5_s_230_n_50: 5252
% t_90_p_5_s_640_n_50: 3457
% -> choose t_90_p_5_s_230_n_50: 5252

n_points = 5252;
n_runs = 100;
grid_size = 101;

max_row_profiles = zeros(n_runs, 2);
min_row_profiles = zeros(n_runs, 2);
max_column_profiles = zeros(n_runs, 2);
min_column_profiles = zeros(n_runs, 2);


for i_runs = 1:n_runs
    % Sample data
    data = copularnd('Gaussian', 0, n_points);
    % Create Empirical CA
    copca_model = CopulaCAEmpirical(data, grid_size);
    % Store profile values
    rowProfiles = copca_model.getRowProfiles(2);
    columnProfiles = copca_model.getColumnProfiles(2);
    for i_dim = 1:2
        max_row_profiles(i_runs, i_dim) = max(rowProfiles(:, i_dim));
        min_row_profiles(i_runs, i_dim) = min(rowProfiles(:, i_dim));
        max_column_profiles(i_runs, i_dim) = max(columnProfiles(:, i_dim));
        min_column_profiles(i_runs, i_dim) = min(columnProfiles(:, i_dim));
    end % i_dim
end % i_runs

save(sprintf('results/bootstrap_independence_profiles/20231026_results_n_%i.mat', grid_size), ...
    'n_points', 'n_runs', 'grid_size', "max_row_profiles", ...
    "min_row_profiles", "max_column_profiles", 'min_column_profiles');
