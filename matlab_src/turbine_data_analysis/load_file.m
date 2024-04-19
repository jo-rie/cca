function [data] = load_file(time, pressure, size, base_path)
%GET_FILENAME get the filename based on parameters from the setup
fname = sprintf('usim_%i%i%i.mat', time, pressure, size);
import_struct = load(fullfile(base_path, fname));
data = import_struct.usim;
end

