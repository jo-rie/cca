function [copula_cdf_funcs] = functionarrayarchimedeancopulas()
    %FUNCTIONARRAYARCHIMEDEANCOPULAS Get a cell array of cdf 3-argument cdf-functions (u, v, theta) for the list of Archimedean copulas from Nelsen book
    S = @(u, v, theta) u + v - 1 - theta * (1./u + 1./v - 1);
    
    
    copula_cdf_funcs = { ...
      @(u, v, theta) max(u.^(-theta) + v.^(-theta) - 1, 0).^(-1 / theta), ... %1
      @(u, v, theta) max(1 - ((1 - u).^theta + (1-v).^theta).^(1/theta), 0), ... %2
      @(u, v, theta) u .* v ./ (1 - theta .* (1 - u) .* (1 - v)), ... %3
      @(u, v, theta) exp(- (((-log(u)).^theta + (- log(v)).^theta).^(1/theta))), ... %4
      @(u, v, theta) - 1 / theta * log(1 + ((exp(-theta * u)- 1).*(exp(-theta * v) - 1) ./ (exp(-theta) - 1))), ... %5
      @(u, v, theta) 1 - ((1-u).^theta + (1-v).^theta - (1-u).^theta.*(1-v).^theta).^(1/theta), ... %6
      @(u, v, theta) max(theta * u .* v + (1-theta) * (u+v-1), 0), ... %7
      @(u, v, theta) max((theta.^2 .* u .* v - (1-u).*(1-v))./(theta.^2 - (theta-1).^2 .*(1-u) .* (1-v)), 0), ... %8
      @(u, v, theta) u .* v .* exp(- theta * log(u) .* log(v)), ... %9
      @(u, v, theta) u .* v ./ (1 + (1-u.^theta) .* (1 - v.^theta)).^(1 / theta), ... %10
      @(u, v, theta) (max(u.^theta .* v.^theta - 2 * (1-u.^theta) .* (1-v.^theta), 0)).^(1/theta), ... %11
      @(u, v, theta) (1 + ((u.^(-1) - 1).^theta + (v.^(-1) - 1).^theta).^(1 / theta)).^(-1), ... %12
      @(u, v, theta) exp(1 - ((1-log(u)).^theta + (1-log(v)).^theta - 1).^(1 / theta)), ... %13
      @(u, v, theta) (1 + ((u.^(-1 / theta) - 1).^theta + (v.^(-1/theta) - 1).^theta).^(1 / theta)).^(-theta), ... %14
      @(u, v, theta) max(1 - ((1-u.^(1/theta)).^theta + (1 - v.^(1/theta)).^theta).^(1/theta), 0).^theta, ... %15
      @(u, v, theta) 0.5 * (S(u, v, theta) + sqrt(S(u, v, theta).^2 + 4 * theta)), ... %16
      @(u, v, theta) (1 + (((1+u).^(-theta) - 1) .* ((1+v).^(-theta) - 1)) ./ (2.^(-theta) - 1)).^(-1/theta) - 1, ... %17
      @(u, v, theta) max(1 + theta ./ log(exp(theta ./(u-1)) + exp(theta ./ (v-1))), 0), ... %18
      @(u, v, theta) theta ./ log(exp(theta./u) + exp(theta ./ v) - exp(theta)), ... %19
    ...   @(u, v, theta) log(exp(u.^(-theta)) + exp(v.^(-theta)) - exp(1)).^(-1 / theta), ... %20
      @(u, v, theta) 1 - (1 - max((1 - (1 - u).^theta).^(1/theta) + (1 - (1-v).^theta).^(1/theta) - 1, 0).^theta).^(1/theta), ... %21
    ...   @(u, v, theta) max((1 - (1-u.^theta) .* sqrt(1 - (1 - v.^theta).^2) - (1 - v.^theta) .* sqrt(1 - (1 - u.^theta).^2)).^(1/theta), 0) ... %22
    };
    end
    
    