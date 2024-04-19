  classdef CopulaMAR < CopulaCA
    %CopulaMAR models additive model following Greenacre using
    
    properties
        eta % parameter eta in model
        discPdfRaw % 
    end
    
    methods
        function obj = CopulaMAR(cdf_func, n, eta)
            %ADDITIVEMODAB Construct an instance of this class
            %   cdf_func: function handle representing the cdf of the
            %   copula
            %   n: grid size
            %   eta: parameter eta in the model
                        
            discPdfRaw = CopulaCAByCdfFunction.evaluate_cdf_func(cdf_func, n);
            
            discMAR = discPdfRaw + ...
                eta * eye(n) - (1 + eta) ./ n .* ones(n, n); % discMAR is centered

            obj@CopulaCA(discMAR, n);
            obj.discPdfRaw = discPdfRaw;
            obj.eta = eta;
        end
        
        function [mat] = getApproxPDFMatrix(obj, k, correct)
            arguments
                obj
                k
                correct = false
            end
            n = obj.n;
            if k > obj.kComputed
                obj = refineK(obj, k);
            end
            mat = obj.u(:, 1:k) * obj.s(1:k, 1:k) * obj.v(:, 1:k)';
            mat = (mat - obj.eta * eye(n) + ...
                (1 + obj.eta) ./ n * ones(n, n));
            if correct
                mat = make_matrix_doubly_stochastic(mat);
            end
        end

    end
end

