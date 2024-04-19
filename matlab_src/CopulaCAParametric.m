classdef CopulaCAParametric < CopulaCAByCdfFunction & AbstractParametricCopulaCA
    %CopulaCAParametric models parametric copulas in terms of their pdf
    
    methods
        function obj = CopulaCAParametric(copFamily, param, n, type)
            arguments
                copFamily
                param
                n
                type = 'Kendall'
            end
            %PCopCAPDF Construct an instance of class PCopCAPDCopulaCAParametricF
            % copFamily: (one-)parametric copula family {'Gaussian', 'Clayton', 'Frank', 'amh', 'Gumbel'}
            % copTau: Kendall's rank correlation
            % n: grid resolution
            
            if lower(type) == "kendall"
                copParam = mycopulaparam(copFamily, param);
                copTau = param;
            elseif lower(type) == "theta"
                copParam = param;
                try
                    copTau = mycopulastat(copFamily, param);
                catch ME
                    warning("Kendall's tau cannot be computet for type " + num2str(type) + ", family " + copFamily + " and param " + num2str(param) + "." + ME.identifier);
                    copTau = NaN;
                end
                    
            else
                error('Type %s not known.', type);
            end
            cdf_func = @(u, v) mycopulacdf(copFamily, ...
                [u, v], copParam);
            
            obj@CopulaCAByCdfFunction(cdf_func, n);
            obj@AbstractParametricCopulaCA(copFamily, copTau);
        end

    end
end

