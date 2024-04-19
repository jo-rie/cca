function [U, d] = eig_sorted(mat)
    % Compute the sorted eigenvalues and eigenvectors of a matrix
    [U, d] = eig(mat);
    d_diag = diag(d);
    [~, ind] = sort(-d_diag);
    d = diag(d_diag(ind));
    U = U(:, ind);
end