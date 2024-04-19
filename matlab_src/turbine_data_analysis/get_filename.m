function [fname] = get_filename(time,pressure,size)
%GET_FILENAME get the filename based on parameters from the setup
fname = sprintf('usim_%i%i%i.mat', time, pressure, size);
end

