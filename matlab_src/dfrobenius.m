function [fro_diff] = dfrobenius(cop_obj1, cop_obj2, k)
    %DFROBENIUS Compute Frobenius distance between the two copula
    %objects; k specifies the maximal degree for which to compute (default is
    %the grid size n)
    arguments
        cop_obj1
        cop_obj2
        k = 0
    end
        if k == 0
            k = cop_obj1.n;
        end
        if k == cop_obj1.n
            fro_diff = norm(cop_obj1.getPDFMatrix() - cop_obj2.getPDFMatrix(), 'fro');
        else
            [U1, S1, V1] = get_usv(cop_obj1, k);
            [U2, S2, V2] = get_usv(cop_obj2, k);
        
            phi1_squared = sum(diag(S1(1:k, 1:k)).^2);
            phi2_squared = sum(diag(S2(1:k, 1:k)).^2);
            phiP_sum = computePhiSquareProductCopula(U1, S1, V1, U2, S2, V2, k);
        
            fro_diff = sqrt(phi1_squared + phi2_squared - 2 * phiP_sum);
        end
    end
    
    
    
    