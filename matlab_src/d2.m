function [fro_dist] = d2(cop_obj1,cop_obj2, k)
    arguments
        cop_obj1
        cop_obj2
        k = NaN
    end
    %D2 Frobenius distance between copulas with normalization by square root of the sum of phi squares
    
    if isnan(k)
        k = cop_obj1.n;
    end
    fro1_squared = sum(cop_obj1.getSingVals(k).^2);
    fro2_squared = sum(cop_obj2.getSingVals(k).^2);
    fro_dist = dfrobenius(cop_obj1, cop_obj2, k) ./ sqrt(fro1_squared + fro2_squared);
    end
    
    