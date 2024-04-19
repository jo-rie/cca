function [phisquare] = computePhiSquareProductCopula(U1, S1, V1, U2, S2, V2, k)
    %computePhiSquareProductCopula Compute phi square of product of centered copulas with svd (U1, S1, V1) and (U2, S2, V2); 
    % k specifies the maximal degree for which to compute
    phisquare = 0;
                
    for i = 1:k
        for j = 1:k
            phisquare = phisquare + S1(i, i) * S2(j, j) * U1(:, i)' * U2(:, j) * V1(:, i)' * V2(:, j);
        end
    end
end