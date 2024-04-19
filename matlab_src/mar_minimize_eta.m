function [eta, fro_er] = mar_minimize_eta(mat, k)
%MAR_MINIMIZE_ETA Find eta such that frobenius norm of resulting matrix
%gets minimal for approximation of degree k; mat should not be centered!
%

% Check that mat is already centered
if sum(mat * ones(length(mat), 1)) < eps
    warning('mat is centered.')
end

res = fminsearch(@(x) mar_fro4eta(x, k, mat), 0.5);
eta = res;
fro_er = mar_fro4eta(eta, k, mat);
end

