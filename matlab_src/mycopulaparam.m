function param = mycopulaparam(family,varargin)
%MYCOPULAPARAM Wrapper around copulaparam that also accepts the 'amh' copulaparam

varargin_raw = varargin;

if nargin > 0
    family = convertStringsToChars(family);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:copulaparam:WrongNumberOfInputs'));
elseif nargin > 2 && strncmpi(varargin{end-1},'type',max(1,numel(varargin{end-1})))
    type = varargin{end};
    types = {'kendall' 'spearman'};
    type = internal.stats.getParamVal(type,types,'''Type''');
    varargin = varargin(1:end-2);
else
    type = 'kendall';
end

families = {'gaussian','t','clayton','frank','gumbel', 'amh', 'ca'};
family = internal.stats.getParamVal(family,families,'FAMILY');

% Already stripped off 'type' args, should only be parameters left.
if isequal(family,'t')
    numParamArgs = 2; % provisionally
else
    numParamArgs = 1;
end
if length(varargin) > numParamArgs
    error(message('stats:copulaparam:UnrecognizedInput'));
end

switch family
case {'gaussian','t','clayton','frank','gumbel'}
    param = copulaparam(family, varargin_raw{:});
case {'amh'}
    if type == "kendall"
        kendalls_tau = varargin{1};
        syms theta;
        tau_theta_equation = 0 == 1 - 2 ./ (3 * theta) - ...
            2 * (1 - theta).^2 ./ (3 * theta.^2) .* log(1 - theta) - kendalls_tau;
        param = vpasolve(tau_theta_equation);
        if isempty(theta) && abs(kendalls_tau) < 10 .* eps
            param = 0; 
        end
    else
        error('Not Implemented');
    end
case {'ca'}
    if type == "kendall"
        kendalls_tau = varargin{1};
        param = 2 * kendalls_tau / (1 + kendalls_tau);
    else
        error('Not Implemented');   
    end
end % switch family


end
