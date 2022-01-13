function assignVars(varargin)
% assignVars('Variable',value,'Var2',value2,...) assigns the values to the variables in the caller function.
if isempty(varargin)
    return;
end
if numel(varargin)==1
    varargin = varargin{1};
end
if mod(numel(varargin),2) ~= 0
    error('There must be an even number of inputs. assignVars(''variable'',value,...)');
end

for i=1:2:numel(varargin)-1,
    if ischar(varargin{i})
        assignin('caller',varargin{i},varargin{i+1});
    else
        error('Invalid variable name. Variable name to assign must be a string')
    end
end
