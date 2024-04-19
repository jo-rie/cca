function p = mycopulacdf(family, u, varargin)
% COPULACDFONEPARAM Wrapper around Matlab's COPULACDF function that accepts more families.
% family {'gaussian', 't', }

if nargin > 0
    family = convertStringsToChars(family);
end

if nargin < 3
    error(message('stats:copulacdf:WrongNumberOfInputs'));
end

[~,d] = size(u);
if d < 2
    error(message('stats:copulacdf:TooFewDimensions'));
end

% Map values outside of unit hypercube to the edges so that they get
% appropriate CDF values.
u(u<0) = 0; % doesn't overwrite NaNs
u(u>1) = 1;

families = {'gaussian', 't', 'clayton', 'frank', 'gumbel', 'amh', 'ca', 'fgm'};
family = internal.stats.getParamVal(family,families,'FAMILY');

% TODO Check parameter ranges

switch family
case {'gaussian','t','clayton','frank','gumbel'}
    p = copulacdf(family, u, varargin{:});
case {'amh'}
    theta = varargin{1};
    if abs(theta) < eps
        p = u(:, 1) .* u(:, 2);
    else
        p = u(:, 1) .* u(:, 2) ./ (1 - theta .* (1 - u(:, 1)) .* (1 - u(:, 2)));
    end
case {'ca'}
    theta = varargin{1};
    p = min([u(:, 1), u(:, 2)], [], 2).^theta .* (u(:, 1) .* u(:, 2)).^(1 - theta);
case {'fgm'}
    theta = varargin{1};
    p = u(:, 1) .* u(:, 2) + theta .* u(:, 1) .* u(:, 2) .* (1 - u(:, 1)) .* (1 - u(:, 2));
end % switch family

end