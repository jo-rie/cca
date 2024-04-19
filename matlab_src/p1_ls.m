function [mat_p1] = p1_ls(mat)
    %% Solve Problem P1 (Frobenius-nearest matrix with row and column sums 1) using the linear system of the KKT equations and thus also for asymmetric matrices
    n = length(mat);
    re_Seite = [2 .* mat(:); ones(2 * n - 1, 1)]; 
    A = zeros(2*n + n.^2 - 1);
    for i = 1:n %rows 1 : n.^2
        for j = 1:n
            row = (i-1) * n + j;
            A(row, row) = 2;
            A(row, n.^2 + j) = 1;
            if i < n
                A(row, n.^2 + n + i) = 1;
            end
        end
    end 
    for i = 1:n % rows n.^2 + 1 : n.^2 + n
        A(n.^2 + i, n * (i-1) + linspace(1, n, n)) = 1;
    end
    for j = 1:n - 1 % rows n.^2 + n + 1 : n.^2 + 2*n
        A(n.^2 + n + j, n * linspace(0, n - 1, n) + j) = 1;
    end
    lgs_solution = A \ re_Seite;
    mat_p1 = reshape(lgs_solution(1:n.^2), [n n]);
end

% Test for SYMMETRIC approximated discrete copula
% n = 5;
% 
% cop_obj = PCopCAPDF('Clayton', .9, n);
% mat = cop_obj.getApproxPDF(3, false);
% mat = mat .* (mat >= 0);
% 
% % Compute Nearest Bistochastic Matrix
% x_numerical = p1_asymmetric_numerical(mat);
% x_ls = p1_ls(mat);
% 
% fprintf('Distance numerical - mat: %.4f\n', norm(mat - x_numerical, 'fro'))
% fprintf('Distance ls - mat: %.4f\n', norm(mat - x_ls, 'fro'))
% fprintf('Distance ls - numerical: %.4f\n', norm(x_ls - x_numerical, 'fro'))
