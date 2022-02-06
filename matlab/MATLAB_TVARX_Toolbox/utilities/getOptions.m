function varargout = getOptions(options,varargin)
% GETOPTIONS: Returns options values in the workspace.
%
% Syntax: 
%
%       [value1,value2,...] = getOptions(options,field1,default1,field2,default2,...)
%
% Description:
%
%       Returns options values in an options structure. Variables with the field names will be created in the caller's 
%       workspace and set to the value in the option variables field (if it exists) or to the default value. If options 
%       contains a field name not in the list passed to getopts, a warning is issued.
%
%       Example:
%           options = struct('tol',1);
%           getOptions(options,'tol',1e-8,'maxits',100);
%           The function have two variable defined in the local workspace, tol with a value of 1 and maxits with a value of 100.
%
% Input Arguments:
%
%       options	-   [struct]    a structure variable
%       field	-   [char]      a field name
%       default	-   [class]     a default value 
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       value	-   [class]     value in the options field (if it exists) or the default value
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       none
%
% References:
%
%       none
%
% Notes:
%
%       none
%
% Copyright:
%
%       (c) Mario J. Miranda & Paul L. Fackler
%           Applied Computational Economics and Finance
%           2002, MIT Press, Cambridge MA.
%           
%           Revision:   1)	Use dynamic fieldnames with structure instead of GETFIELD, by Gabriel Bruneau, 2013


% Input and Output arguments checking
% ___________________________________

narginchk(1,inf); nargoutchk(0,1);

validateattributes(options,{'struct'},{},'getOptions','options',1)

if isstruct(options)
	optstruct = 1;
else
	optstruct = 0;
end

K = fix(nargin/2);
if nargin/2 == K
	error('getOptions:PairingFieldsValues','getOptions: Fields and default values must come in pairs')
end


% Get options
% ___________

varargout = cell(K,1);
k = 0;

ii = 1;
for i = 1:K
	if optstruct && isfield(options,varargin{ii})
        assignin('caller',varargin{ii},options.(varargin{ii})); % Use dynamic fieldnames with structure instead of GETFIELD
        k = k+1;
	else
        assignin('caller',varargin{ii},varargin{ii+1});
	end
	ii = ii + 2;
end
  
if optstruct && k ~= size(fieldnames(options),1)
	warning('getOptions:ImproperFields','getOptions: Options variables contains improper fields')
end


end % function {getOptions}

