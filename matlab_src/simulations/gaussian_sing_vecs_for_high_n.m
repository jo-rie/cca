function gaussian_sing_vecs_for_high_n(ns, params, degrees, exp_path)
%GAUSSIAN_SING_VECS_FOR_HIGH_N Compute the singular vectors for Gaussian
%copula and high ns
for i_ns = 1:length(ns)
    n = ns(i_ns);
    fprintf('Starting n = %i\n', n);
    for i_tau = 1:length(params)
        singular_vectors = zeros([n, degrees]);
        param = params(i_tau);
        cop_obj = CopulaCAParametric('Gaussian', param, n, 'theta');
        singular_vectors(:, :) = cop_obj.getU(degrees);
        save( ...
          fullfile( ...
            exp_path, ...
            sprintf( ...
              'gaussian_sing_vecs_for_n_%i_rho_%.2f.mat', n, param)), ...
          'n', 'param', 'degrees', 'singular_vectors' ...
        )
    end
    
end
end

