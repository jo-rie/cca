classdef CopulaCA<handle
    % CopulaCA models a discrete copula through a singular value
    % decomposition
    
    properties (Access = public)
        n {mustBePositive} = 1
        plotGrid
        
        discCenteredPDF double % discretized PDF matrix
        discCDF double % discretized CDF matrix
        
        u % matrix of left singular values
        v % matrix of right singular values
        s % matrix with singular values on the diagonal

        kComputed int16 % number of singular values computed
    end

    properties (Access = public, Constant = true)
        TOL = 1e-7
    end
    
    methods
        function obj = CopulaCA(discCenteredPDF, n)
            % COPULACA Constructor of class
            obj.n = n;
            obj.discCenteredPDF = discCenteredPDF;

            obj.plotGrid = linspace(1 / n / 2, 1 - 1 / n / 2, n);      
            obj.discCDF = CopulaCA.pdf2cdf(obj.discCenteredPDF);

            obj.kComputed = min([obj.n, 20]);

            obj = refineK(obj, obj.kComputed);
        end
        
        % TODO add function to evaluate pdf and cdf at a point (u, v)

        function mat = getPDFMatrix(obj)
            % Returns the full pdf matrix
            mat = obj.discCenteredPDF + ones(obj.n, obj.n) ./ obj.n;
        end
        
        function mat = getCDFMatrix(obj)
            % Return the full cdf matrix
            mat = obj.discCDF;
        end
        
        function mat = getU(obj, k)
            % Matrix of left singular values up to rank k
            if k <= obj.kComputed
                mat = obj.u(:, 1:k);
            else
                refineK(obj, k);
                mat = obj.u(:, 1:k);
            end
        end
        
        function mat = getV(obj, k)
            % Matrix of right singular values up to rank k
            if k <= obj.kComputed
                mat = obj.v(:, 1:k);
            else
                refineK(obj, k);
                mat = obj.v(:, 1:k);
            end
        end
        
        function mat = getRowProfiles(obj, k)
            % Returns the correspondence analysis row profiles up to rank k
            arguments
                obj
                k = 2
            end
            mat = obj.u(:, 1:k) * obj.s(1:k, 1:k); 
        end

        function mat = getColumnProfiles(obj, k)
            % Returns the correspondence analysis column profiles up to rank k
            arguments
                obj
                k = 2
            end
            mat = obj.v(:, 1:k) * obj.s(1:k, 1:k); 
        end

        function mat = getS(obj, k)
            % Returns a n-times-n matrix with the largest k singular values on the diagonal
            if k <= obj.kComputed
                mat = obj.s(1:k, 1:k); 
            else
                refineK(obj, k);
                mat = obj.s(1:k, 1:k); 
            end
        end

        function singvals = getSingVals(obj, k)
            % Returns a vector with the largest k singular values on the diagonal
            singvals = diag(obj.getS(k));
        end

        function [res] = getPDFMatrixApproxError(obj, k, correct)
            arguments
                obj
                k
                correct = false
            end
            % Compute the approximation error of the PDF matrix when using a rank-k-approximation. correct specifies whether the approximation should be corrected for negative elements.
            % The result is a structure with different matrix norm results.
            res = struct();
            pdfApprox = obj.getApproxPDFMatrix(k, correct);
            pdfError = obj.getPDFMatrix() - pdfApprox;
            res.('l1') = norm(pdfError, 1);
            res.('l2') = norm(pdfError, 2);
            res.('linf') = norm(pdfError, Inf);
            res.('fro') = norm(pdfError, 'fro');
        end

        function [res] = getCenteredPDFMatrix(obj)
            res = obj.discCenteredPDF;
        end

        function [res] = getCDFApproxError(obj, k, correct)
            arguments
                obj
                k
                correct = false
            end
            % Compute the approximation error of the CDF matrix when using a rank-k-approximation. correct specifies whether the approximation should be corrected for negative elements.
            % The result is a structure with different matrix norm results.
            res = struct();
            cdfApprox = obj.getApproxCDFMatrix(k, correct);
            cdfError = obj.discCDF - cdfApprox;
            res.('l1') = norm(cdfError, 1);
            res.('l2') = norm(cdfError, 2);
            res.('linf') = norm(cdfError, Inf);
            res.('fro') = norm(cdfError, 'fro');
        end

        function [mat, meaning] = getInvalidCDFValues(obj, k)
            % Get a matrix of indicating invalid CDF values in the rank k approximation. meaning specifies the encoding of the matrix
            meaning=["valid", "row/col sum != u/v", "neg_value", "not 2-inc"];
            mat = zeros(obj.n, obj.n);

            approx = getApproxCDFMatrix(obj, k, false);
            
            mat(abs(sum(approx, 1) - 1/obj.n:1/obj.n:1) > obj.TOL, end) = 1;
            mat(end, abs(sum(approx, 2) - 1/obj.n:1/obj.n:1) > obj.TOL) = 1;
            
            mat(approx < 0) = 2;
            
            mat(1:end-1, 1:end-1) = (approx(2:end, 2:end)...
                - approx(1:end-1, 2:end) - approx(2:end, 1:end-1)...
                + approx(1:end-1, 1:end-1) < 0) * 3;
        end

        function [mat, meaning] = getInvalidPDFValues(obj, k)
            % Get a matrix of indicating invalid PDF values in the rank k approximation. meaning specifies the encoding of the matrix
            meaning=["valid", "row/col sum != 1", "neg. entry"];
            mat = zeros(obj.n, obj.n);
            
            approx = getApproxPDFMatrix(obj, k, false);
            
            mat(abs(sum(approx, 1) - 1) > 1e-5, end) = 1;
            mat(end, abs(sum(approx, 2) - 1) > 1e-5) = 1;
            
            mat(approx < 0) = 2;
        end

        function [geom_dim] = getGeometricDimension(obj, tol)
            % Compute the geometric dimension of the object.
            arguments
                obj
                tol = 1e-15
            end
            if isnan(tol)
                geom_dim = rank(obj.discCenteredPDF);
            else
                geom_dim = rank(obj.discCenteredPDF, tol);
            end
        end

        function [tau] = getKendallsTau(obj, k)
            % Compute Kendalls tau of an approximation of rank k of the object.
            arguments
                obj
                k = obj.n
            end
            % Use Matrix formulation
            E = tril(2 * ones(obj.n)) - eye(obj.n);
            tau = 1 - 1 ./ obj.n .^2 * ...
                trace(E * obj.getApproxPDFMatrix(k, false) * ...
                E * obj.getApproxPDFMatrix(k, false)');
        end

        function [rho] = getSpearmansRho(obj, k, correct)
            % Compute Spearmans rho of an approximation of rank k of the object.
            arguments
                obj
                k = obj.n
                correct = false
            end
            % Use Matrix formulation
            omega = (2 .* obj.n + 1 - 2 * (1:obj.n)) ./ obj.n;
            Omega = omega' * omega;
            rho = 3 ./obj.n .* trace(Omega * ...
                obj.getApproxPDFMatrix(k, correct)) - 3;
        end

        function [obj] = refineK(obj, k)
            % Compute a finer SVD of order k
            if k > obj.n
                warning('Attempting to compute a finder SVD than grid size. Aborting')
            else
                obj.kComputed = k;
                [obj.u, obj.s, obj.v] = svds(obj.discCenteredPDF, k);
                % Chose singular vectors such that smallest non-zero element is
                % positive
                % Get first non-zero element by X(find(X, 1))
                for j = 1:size(obj.u, 2)
                    if obj.u(find(obj.u(:, j), 1), j) < 0
                        obj.u(:, j) = obj.u(:, j) * (-1);
    %                 end % Remove, as the approximation changes if only one
    %                 singular vector is multiplied by 1
    %                 if obj.v(find(obj.v(:, j), 1), j) < 0
                        obj.v(:, j) = obj.v(:, j) * (-1);
                    end
                end
            end
        end

        function [mat] = getApproxPDFMatrix(obj, k, correct)
            arguments
                obj
                k
                correct = false
            end
            % Get approximated PDF of order k. If correct, the approximation's negative elements are corrected by computing the (Frobenius-)nearest doubly stochastic.
            if k > obj.kComputed
                obj = refineK(obj, k);
            end
            mat = obj.u(:, 1:k) * obj.s(1:k, 1:k) * obj.v(:, 1:k)';
            mat = mat + ones(obj.n, obj.n) / obj.n;
            if correct
                mat = make_matrix_doubly_stochastic(mat);
            end
        end

        function [mat] = getApproxCDFMatrix(obj, k, correct)
            arguments
                obj
                k
                correct = false
            end
            pdfmat = getApproxCDF(obj, k, correct);
            mat = pdf2cdf(pdfmat);
        end

        function [u, s, v] = get_usv(obj, k)
            if k > obj.kComputed
                refineK(obj, k);
            end
            u = obj.u;
            s = obj.s;
            v = obj.v;
        end
    end
    
    methods (Static = true)
        function [mat] = pdf2cdf(pdfMat)
            % Compute the cdf matrix based on the given pdf matrix
            mat = zeros(size(pdfMat) + 1);
            mat(2:end, 2:end) = cumsum(cumsum(pdfMat)');
        end
        
        function [mat] = cdf2pdf(cdfMat)
            % Compute the pfd matrix based on the given cdf matrix
            mat = (cdfMat(2:end, 2:end) - cdfMat(1:end-1, 2:end) ...
                - cdfMat(2:end, 1:end-1) + cdfMat(1:end-1, 1:end-1)) .* length(cdfMat);
            % Problem: pdf matrix has size n - 1
        end
    end
    
end

