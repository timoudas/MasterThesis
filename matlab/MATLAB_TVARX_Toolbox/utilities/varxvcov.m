function vcov = varxvcov(x,resid,het,corre)
% VARXVCOV: Estimate the parameters' variance-covariance matrix estimated using VARX.
%
% Syntax:
%
%       vcov = varxvcov(x,errors,het,uncorr)
%
% Description:
%
%       Estimate the variance-covariance matrix for parameters estimated using VARX and TVARX under one of these
%       conditions 
%           1 - Conditionally Homoskedastic and Uncorrelated Errors
%           2 - Conditionally Homoskedastic but Correlated Errors
%           3 - Heteroskedastic but Conditionally Uncorrelated Errors
%           4 - Heteroskedastic and Correlated Errors
%
% Input Arguments:
%
%       x       -	[double]    a numObs-by-numX matrix of regressors
%       resid   -   [double]    a numObs-by-numY vector of residuals
%       het     -   [logical]   type of covariance estimator
%                                   true    heteroskedastic
%                                   false	homoskedastic
%       corr	-   [logical]   the assumed structure of the error covariance matrix
%                                   true    correlated errors
%                                   false	uncorrelated errors
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       vcov	-   [double]    a (numY*numX)-by-(numY*numX) matrix of estimated parameter variance-covariance
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
%       Helper function for VARX and TVARX. 
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014
%           based (heavily) on vectorarvcv.m by
%               Kevin Sheppard (kevin.sheppard@economics.ox.ac.uk)


% Input and Output arguments checking
% ___________________________________

narginchk(4,4); nargoutchk(1,1);

[x,resid,het,corre] = checkInput(x,resid,het,corre);


% Define sizes and constants
% __________________________

[numObs,numVars] = size(resid);
numX = size(x,2);

sige = resid'*resid/numObs;
xpxi = ((x'*x)/numObs)^(-1);
x2 = repmat(x,1,numVars);
e2 = reshape(repmat(resid,numX,1),numObs,numVars*numX);
s = x2.*e2;


% Variance-covariance estimation
% ______________________________

switch het
    case 1  % Heteroskedastic residuals
        switch corre 
            case 1  % Correlated residuals
                % White variance-covariance matrix
                Ainv = kron(eye(numVars),xpxi);
                B = (s'*s)/numObs;
                vcov = Ainv*B*Ainv/numObs;                
            
            case 0	% Uncorrelated residuals
                Ainv = kron(eye(numVars),xpxi);    
                B = zeros(numX*numVars);
                for i = 1:numVars
                    sel = (i-1)*numX+1:i*numX;
                    temp = s(:,sel);
                    B(sel,sel) = temp'*temp/numObs;
                end
                vcov = Ainv*B*Ainv/numObs;        
        
        end % switch {corre}
    
    case 0	% Homoskedastic residuals
        switch corre 
            case 1  % Correlated residuals
                vcov = kron(sige,xpxi)/numObs;                
            
            case 0	% Uncorrelated residuals
                vcov = kron(diag(diag(sige)),xpxi)/numObs;        
        
        end % switch {corre}
end % switch {het}
          

end % function {varxvcov}


% ---------------------------------
function [x,resid,het,corre] = checkInput(x,resid,het,corre)
% checkInput: Local function to check the validity of required inputs

% x
validateattributes(x,{'numeric'},{'real','2d','real','finite','nonnan'},'varxvcov','x',1);

% resid
validateattributes(resid,{'numeric'},{'real','2d','real','finite','nonnan'},'varxvcov','resid',2);

% het
validateattributes(het,{'logical'},{'numel',1},'varxvcov','het',3);

% corre
validateattributes(corre,{'logical'},{'numel',1},'varxvcov','corre',4);


end % subfunction {checkInput}

