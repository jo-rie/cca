function counterexample_polygon_runs(exp_path, center_points, polyshapes, ns)
    cdf_func = @(u, v) copulacdf_polygons(u, v, polyshapes);
    
    geom_dims = zeros(length(ns), 1);

    for i_n = 1:length(ns)
        cop_obj = CopulaCAByCdfFunction(cdf_func, ns(i_n));
        geom_dims(i_n) = cop_obj.getGeometricDimension(1e-8);
    end %i_n

    save(fullfile(exp_path, 'geom_dims.mat'), 'ns', 'center_points', 'geom_dims')
end %counterexample_polygon_runs