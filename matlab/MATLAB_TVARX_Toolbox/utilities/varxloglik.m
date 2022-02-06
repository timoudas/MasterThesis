function [logL,logCL] = varxloglik(resid)
% VARXLOGLIK: Computes total and conditional loglikelihoods of a multidimensional VARX time series process.
%
% Syntax:
%
%       [logL,logCL] = varxloglik(resid)
%
% Description:
%
%       VARXLOGLIK computes total and conditional loglikelihoods of a multidimensional VARX 
%       time series process.
%
% Input Arguments:
%
%       resid       -   [double]	(numObs-numNaN)-by-numY, residuals
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       logL        -	[double]	scalar, the total loglikelihood of the
%                                   response data in each path of W. logL = sum(logCL).
%
%       logCL       -	[double]    numObs-by-1 matrix containing the conditional
%                                   loglikelihoods of the response data based on residuals.
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
%           based on vgxloglik.m of the MATLAB Econometrics Toolbox


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(1,2);

resid = checkInput(resid);


% Define constant
% _______________

numObs = size(resid,1);
numEqs = size(resid,2);


% Residuals variance-covariance matrix
% ____________________________________

sige = resid.'*resid/numObs;


% Compute the log-likelihood
% __________________________

% Conditional log-likelihood
logCL = zeros(numObs,1);
LLP = log(det(sige));
for t = 1:numObs
	LLQ = resid(t,:)/sige*resid(t,:).'; 
	logCL(t) = -0.5*(numEqs*log(2*pi) + LLP + LLQ);
end

% Total log-likelihood
logL = -(numEqs*numObs/2)*log(2*pi) - numObs/2 * log(det(sige)) - 1/2*resid(:).'*kron(inv(sige),eye(numObs))*resid(:);


end % function {varxloglik}


% ---------------------------------
function resid = checkInput(resid)
% checkInput: Local function to check the validity of required inputs

% resid
validateattributes(resid,{'numeric'},{'real','2d','finite','nonnan'},'varxloglik','resid',1);  % check resid


end % subfunction {checkInput}



