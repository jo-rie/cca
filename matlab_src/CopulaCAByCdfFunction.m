classdef CopulaCAByCdfFunction < CopulaCA
    %CopulaCAByCdfFunction models a copula given by its CDF function
    
    methods
        function obj = CopulaCAByCdfFunction(cdf_func, n)
            %PCopCAPDF Construct an instance of class CopulaCAByCdfFunction
            % cdf_func: copula cdf (function taking two arguments (u, v))
            % n: grid size

            % Create frequency matrix
            discPDF = CopulaCAByCdfFunction.evaluate_cdf_func(cdf_func, n);
            discCenteredPDF = discPDF - ones(n, n) ./ n;
            obj@CopulaCA(discCenteredPDF, n);
        end

    end

    methods (Static)
        function discPDF = evaluate_cdf_func(cdf_func, n)
            % Compute the discretized copula PDF matrix for the given cdf_func with grid size n.
            % cdf_func should be a two-argument CDF function that accepts two column vectors.
            grid = linspace(eps, 1-eps, n+1); % Compute the discrete cdf matrix for grid size n+1
            [gridx, gridy] = meshgrid(grid);
            values = cdf_func(gridx(:), gridy(:));
            values = reshape(values, n + 1, n + 1);
            discPDF = (values(2:end, 2:end) - values(1:end-1, 2:end) ...
                - values(2:end, 1:end-1) + values(1:end-1, 1:end-1)) .* n;
        end
    end
end

