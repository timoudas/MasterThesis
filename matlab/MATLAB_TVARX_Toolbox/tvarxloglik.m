function [logL,logCL] = tvarxloglik(resid,commonvar)
% TVARXLOGLIK: Computes total and conditional loglikelihoods of a multidimensional TVARX time series process.
%
% Syntax:
%
%       [logL,logCL] = tvarxloglik(resid,commonvar)
%
% Description:
%
%       TVARXLOGLIK computes total and conditional loglikelihoods of a multidimensional TVARX time series process.
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
%       logL        -	[double]	scalar, the total loglikelihood of the response data in each path of W. logL = sum(logCL).
%       logCL       -	[double]    numRegimes-by-1 cell array, with each cell containing a numObs-by-1 matrix containing the conditional
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


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(1,2);

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
        % Compute the conditional log-likelihood
        % ______________________________________

        logCL = cell(size(resid,1),1);

        for indRegimes = 1:size(resid,1)

            % Regimes residuals
            residReg = resid{indRegimes};
            validateattributes(residReg,{'numeric'},{'real','2d','finite'},'tvarxloglik','resid{indRegimes}',1);  % check resid{indRegimes}

            numObsReg = size(residReg,1);

            % Residuals variance-covariance matrix
            sige = residReg.'*residReg/numObsReg;

            % Log-likelihood
            logCL{indRegimes} = zeros(numObsReg,1);

            LLC = size(residReg,2)*log(2*pi);
            LLP = log(det(sige));
            for t = 1:numObsReg
                LLQ = residReg(t,:)/sige*residReg(t,:).'; 
                logCL{indRegimes}(t) = -0.5*(LLC + LLP + LLQ);
            end

        end % for {indRegimes}

        % Concatenate all log-likelihood for all regimes
        logLtemp = [];
        for indRegimes = 1:size(resid,1)
            logLtemp = [logLtemp;logCL{indRegimes}]; %#ok<AGROW>
        end % for {indRegimes}


        % Compute the log-likelihood
        % __________________________

        logL = sum(logLtemp);

    case 1
        
        % Residuals variance-covariance matrix
        % ____________________________________
        
        residReg = [];
        for indRegimes = 1:size(resid,1)
            residReg = [residReg;resid{indRegimes}]; %#ok<AGROW>
        end
        
        numObsReg = size(residReg,1);
        sige = residReg.'*residReg/numObsReg;
        
        % Compute the conditional log-likelihood
        % ______________________________________

        logCL = cell(size(resid,1),1);

        for indRegimes = 1:size(resid,1)

            % Regimes residuals
            residReg = resid{indRegimes};
            validateattributes(residReg,{'numeric'},{'real','2d','finite'},'tvarxloglik','resid{indRegimes}',1);  % check resid{indRegimes}

            numObsReg = size(residReg,1);

            % Log-likelihood
            logCL{indRegimes} = zeros(numObsReg,1);

            LLC = size(residReg,2)*log(2*pi);
            LLP = log(det(sige));
            for t = 1:numObsReg
                LLQ = residReg(t,:)/sige*residReg(t,:).'; 
                logCL{indRegimes}(t) = -0.5*(LLC + LLP + LLQ);
            end

        end % for {indRegimes}

        % Concatenate all log-likelihood for all regimes
        logLtemp = [];
        for indRegimes = 1:size(resid,1)
            logLtemp = [logLtemp;logCL{indRegimes}]; %#ok<AGROW>
        end % for {indRegimes}


        % Compute the log-likelihood
        % __________________________

        logL = sum(logLtemp);

        
end % switch {commonvar}


end % function {tvarxloglik}


% ---------------------------------
function resid = checkInput(resid)
% checkInput: Local function to check the validity of required inputs

% resid
validateattributes(resid,{'cell'},{},'tvarxlogdet','resid',1);  % check resid
for indRegimes = 1:numel(resid)
    validateattributes(resid{indRegimes},{'numeric'},{'real','2d','finite','nonnan'},'tvarxlogdet','resid',1);  % check resid
end


end % subfunction {checkInput}

