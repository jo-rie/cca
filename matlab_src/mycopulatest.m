

families = {'gaussian', 'gaussian', 't', 'clayton', 'frank', 'gumbel', 'amh', 'ca', 'fgm'};
theta_families = {0.9, 0.1, [.9, .5], 10, .9, 20, .3, 0.9, .3};

figure;
tiledlayout(2, 5);
for i_copula = 1:length(families)
    nexttile;
    theta = theta_families{i_copula};
    family = families{i_copula};
    u1 = linspace(1e-3,1-1e-3,50);
    [U1,U2] = meshgrid(u1,u1);
    if families{i_copula} == "t"
        f = mycopulacdf(family, [U1(:) U2(:)], theta(1), theta(2));
    else
        f = mycopulacdf(family, [U1(:) U2(:)], theta);
    end
    f = reshape(f,size(U1));
    surf(u1,u1,f,'FaceColor','interp','EdgeColor','none')
    xlabel('u')
    ylabel('v')
    zlabel('CDF')
    if family == "t"
        title(sprintf('%s, %.2f', family, theta(1)))
    else
        title(sprintf('%s, %.2f', family, theta))
    end
end %i_copula
    