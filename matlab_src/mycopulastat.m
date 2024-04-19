function r = mycopulastat(family,varargin)
% MYCOPULASTAT Wrapper around Matlab's copulastat that also accepts the amh and CA copula

varargin_raw = varargin;

if nargin > 0
    family = convertStringsToChars(family);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:copulastat:WrongNumberOfInputs'));
elseif nargin > 2 && strncmpi(varargin{end-1},'type',max(1,numel(varargin{end-1})))
    type = varargin{end};
    types = {'kendall' 'spearman'};
    type = internal.stats.getParamVal(type,types,'''Type''');
    varargin = varargin(1:end-2);
else
    type = 'kendall';
end

families = {'gaussian','t','clayton','frank','gumbel', 'amh', 'ca', 'fgm'};
family = internal.stats.getParamVal(family,families,'FAMILY');

% Already stripped off 'type' args, should only be parameters left.
if isequal(family,'t')
    numParamArgs = 2; % provisionally
else
    numParamArgs = 1;
end
if length(varargin) > numParamArgs
    error(message('stats:copulastat:UnrecognizedInput'));
end

switch family
case {'gaussian','t','clayton','frank','gumbel'}
    r = copulastat(family, varargin_raw{:});
case {'amh'}
    theta = varargin{1};
    if type == "kendall"
        r = 1 - 2 ./ (3 * theta) - 2 * (1 - theta).^2 ./ (3 * theta.^2) .* log(1 - theta);
    else
        error(message('stats:copulastat:TypeNotImplemented'));
    end
case {'ca'}
    theta = varargin{1};
    if type == "kendall"
        r = theta / (2 - theta);
    else
        error(message('stats:copulastat:TypeNotImplemented'));
    end        
case {'fgm'}
    r = nan;
end % switch family

end
