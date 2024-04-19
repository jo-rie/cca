function [fro_dist] = d1(cop_obj1,cop_obj2, k)
    arguments
        cop_obj1
        cop_obj2
        k = NaN
    end
    %D1 Frobenius distance between copulas with normalization by sqrt(2n - 2)
    
    if isnan(k)
        k = cop_obj1.n;
    end
    n = cop_obj1.n;
    fro_dist = dfrobenius(cop_obj1, cop_obj2, k) ./ sqrt(2 * n - 2);
    end
    
    