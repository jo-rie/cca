function [mat_p1] = p1_numerical(mat)
    % Solve problem p1 numericaly using Matlab's fmincon function

    n = length(mat);
    x0 = eye(n, n);
    x0 = x0(:);
    A = zeros(2*n, n.^2);

    for i = 1:n
        A_tmp = zeros(n,n);
        A_tmp(i, :) = 1;
        A_tmp_trapo = A_tmp';
        A(i, :) = A_tmp(:);
        A(n + i, :) = A_tmp_trapo(:);
    end
    
    b = ones(2*n, 1);
    
    res = fmincon(@(x) sum((x - mat(:)).^2), x0, [], [], A, b, [], [], ...
        [], optimoptions('fmincon', 'Display', 'off'));
    mat_p1 = reshape(res, n, []);

    % disp(A * mat_p1(:))
end