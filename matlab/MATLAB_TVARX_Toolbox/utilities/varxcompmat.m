function CompMat = varxcompmat(results)
% VARXCOMPMAT: Construct VAR(1) companion form of a VARX(p) model.
%
% Syntax:
%
%       CompMat = varxcompmat(results)
%
% Description:
%
%       Construct VAR(1) companion form of a VARX(p) model.
%
% Input Arguments:
%
%       results     -   [struct]	a results structure returned by VARX function
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       CompMat     -   [double]    the companion matrix
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput
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
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,1);

results = checkInput(results);


% Compute the companion form matrix
% _________________________________

numY    = results.sizes.numY;
numLags = results.sizes.numLags;

CompMat = zeros(numY*numLags);
for i = 1:numLags
    CompMat(1:numY,numY*i-numY+1:numY*i) = results.beta.lags{i};
end
CompMat(numY+1:end,1:end) = [eye(numY*(numLags-1)) zeros(numY*(numLags-1),numY)];


end % function {varxcompmat}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'varxcompmat','results',1);
if ~strcmpi(results.meth,'varx')
    error('varxcompmat:InputError','varxcompmat: The results structure provided must be returned by VARX function')
end


end % subfunction {checkInput}

