classdef StudentTCopulaCA < CopulaCA
    %StudentTCopulaCA models parametric t copulas in terms of their pdf
    properties
        rho
        nu
    end

    methods
        function obj = TCopCAPDF(rho, nu, n)
            %TCopCStudentTCopulaCA Construct an instance of class TCopCAPDF

            cdf_func = @(u, v) copulacdf('t', ...
                    [u, v], rho, nu);

            
            copFamily = 't';
            copTau = copulastat('t', rho, nu);

            % Call super contructors
            obj@FuncCopCAPDF(cdf_func, n);
            obj@PCopCA(copFamily, copTau);

            % Store T copula parameters
            obj.rho = rho;
            obj.nu = nu;
        end

        function rho = getRho(obj)
            rho = obj.rho;
        end

        function nu = getNu(obj)
            nu = obj.nu;
        end

    end
end

