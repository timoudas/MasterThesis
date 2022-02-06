function logDet = tvarxlogdet(resid,commonvar)
% TVARXLOGDET: Computes log determinant of a multidimensional TVARX time series process.
%
% Syntax:
%
%       logDet = tvarxlogdet(resid,commonvar)
%
% Description:
%
%       TVARXLOGDET computes log determinant of a multidimensional TVARX time series process.
%
% Input Arguments:
%
%       resid       -   [cell]      numRegimes-by-1 cell array, with each cell containing a (numObs-numNaN)-by-numY matrix of residuals
%       commonvar	-   [logical]	each regime share the same distribution
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       logDet      -	[double]	scalar, the log determinant
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

narginchk(1,2); nargoutchk(1,1);

resid = checkInput(resid);


% Options
% _______

if nargin == 2 
    validateattributes(commonvar,{'logical'},{'numel',1},'tvarxlogdet','commonvar',2);
else
    commonvar = false;
end


% Log-likelihood
% ______________

switch commonvar
    case 0
        logDet = 0;
        for indRegimes = 1:size(resid,1)
            sige = (resid{indRegimes}.'*resid{indRegimes})./size(resid{indRegimes},1);
            logDet = logDet + size(resid{indRegimes},1)/2 * log(det(sige));
        end
        
	case 1
        residall = [];
        for indRegimes = 1:size(resid,1)
        	residall = [residall;resid{indRegimes}]; %#ok<AGROW>
        end
        sige = (residall.'*residall)./size(residall,1);
        logDet = log(det(sige));
        
end


end % function {tvarxlogdet}


% ---------------------------------
function resid = checkInput(resid)
% checkInput: Local function to check the validity of required inputs

% resid
validateattributes(resid,{'cell'},{},'tvarxlogdet','resid',1);  % check resid
for indRegimes = 1:numel(resid)
    validateattributes(resid{indRegimes},{'numeric'},{'real','2d','finite','nonnan'},'tvarxlogdet','resid',1);  % check resid
end


end % subfunction {checkInput}

