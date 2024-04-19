function [cdf_vals] = copulacdf_polygons(u, v, polygons)
%COPULACDF_POLYGONS Returns the cdf values at the points (u, v) (u and v
%columns vectors) for a copula with uniform mass on the polygons specified
%as polygon matrix 'polygons'
total_area = 0;
    for k = 1:length(polygons)
        total_area = total_area + area(polygons(k));
    end
    cdf_vals = zeros(length(u), 1);
    
    for i = 1:length(u)
        area_tmp = 0;
        if u(i) > eps && v(i) > eps
            for k = 1:length(polygons)
                area_tmp = area_tmp + area(intersect(...
                    polyshape([0 0 u(i) u(i)], [0 v(i) v(i) 0]), polygons(k)));
            end
        end
        cdf_vals(i) = area_tmp / total_area;
    end
end

