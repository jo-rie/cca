function [doubly_stochastic_mat] = make_matrix_doubly_stochastic(mat, tol, maxIter)
% Compute the frobenius-nearest doubly stochastic matrix to mat based on the algorithm by Zass (2007) and the asymmetric extension. 
% tol specifies the tolerance for with an element >= - tol should be considered zero. 
% maxIter specifies the maximal number of iterations.
arguments
    mat
    tol = 1e-7
    maxIter = 1e5
end



n = length(mat);
if min(mat, [], 'all') >= - tol
    doubly_stochastic_mat = mat;
else 
    % Solve P1 explicitely or by the solution of a linear system
    if all(mat == mat')
        p1 = @(x) x + ...
                (1 / n * eye(n) ...
                + ones(1, n) * x * ones(n, 1) * eye(n, n) ./ n.^2 ... 
                - x / n) * ones(n, n) ...
                - ones(n, n) * x / n;
    else 
        p1 = @(x) p1_ls(x);
    end

    doubly_stochastic_mat = p1(mat);  
    counter = 0;
    while any(doubly_stochastic_mat < - tol, 'all')
        if counter > maxIter
            warning("End of frobenius-optimization procedure after 1e5 runs")
            break
        end
        doubly_stochastic_mat = p2(doubly_stochastic_mat); 
        doubly_stochastic_mat = p1(doubly_stochastic_mat);
        counter = counter + 1;
    end
end

end % function

function [b] = p2(mat)
    % Solution for problem P2
    b = (mat >= 0) .* mat;
end