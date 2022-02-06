function CompMat = tvarxcompmat(results)
% TVARXCOMPMAT: Construct TVAR(1) companion form of a TVARX(p) model.
%
% Syntax:
%
%       CompMat = tvarxcompmat(results)
%
% Description:
%
%       Construct TVAR(1) companion form of a TVARX(p) model.
%
% Input Arguments:
%
%       results     -   [struct]	a results structure returned by TVARX function
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


% Compute the companion form matrices
% _________________________________

CompMat = cell(size(results,1),1);
for indRegimes = 1:size(results,1)
    CompMat{indRegimes} = zeros(results(indRegimes).sizes.numY*results(indRegimes).sizes.numLags);
    for i = 1:results(indRegimes).sizes.numLags
        CompMat{indRegimes}(1:results(indRegimes).sizes.numY,results(indRegimes).sizes.numY*i-results(indRegimes).sizes.numY+1:results(indRegimes).sizes.numY*i) = results(indRegimes).beta.lags{i};
    end
    CompMat{indRegimes}(results(indRegimes).sizes.numY+1:end,1:end) = ...
        [eye(results(indRegimes).sizes.numY*(results(indRegimes).sizes.numLags-1)) zeros(results(indRegimes).sizes.numY*(results(indRegimes).sizes.numLags-1),results(indRegimes).sizes.numY)];
end


end % function {tvarxcompmat}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxcompmat','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxcompmat:InputError','tvarxcompmat: The results structure provided must be returned by TVARX function')
end


end % subfunction {checkInput}

