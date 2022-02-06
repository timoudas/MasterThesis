function infocrit = infocriterion(mode,data,numParams,numObs)
% INFOCRITERION: Compute information criteria for model selection.
% 
% Syntax:
%
%       infocrit = infocriterion(mode,data,numParams,numObs)
%
% Description:
%
%       Compute Akaike (AIC), Corrected Akaike (AICc), Akaike-Bayesian (ABIC), Schwarz-Bayesian (SBIC), Hannan-Quinn (HQIC) 
%       and Akaike Final Prediction Error (FPE) information criteria for order (or parsimony) selection.
%
%       Given optimized log-likelihood function values logL obtained by fitting a model to data or a matrix of residuals 
%       from regression, compute the information criteria. Since information criteria penalize models with additional 
%       parameters, information criteria select models based on both goodness of fit and parsimony. When using either 
%       information criteria, models that minimize the criteria are preferred.
%
%       If numObs/numParams <= 40, then we suggest to use AICc rather then AIC.
%
% Input Arguments:
%
%       mode            -   [char]      string, characterization of the input given to compute the criterion:
%                                           'logL'      for log-likelihood (univariate or multivariate models)
%                                           'ResUni'    for residuals from univariate models
%                                           'ResMulti'	for residuals from multivariate models
%       data            -   [double]    array, with inputs that depend on the "mode":
%                                           If mode = 'logL', a numModels by 1 vector of optimized log-likelihood objective 
%                                               function values associated with parameter estimates of various models.
%                                           If mode = 'ResUni', a numModels by numObs matrix of residuals from various 
%                                               univariate models.
%                                           If mode = 'ResMulti', a numModels by numObs by numVars array of residuals from a 
%                                               various multivariate models.
%       numParams       -   [integer]	vector, numModels by 1 vector, number of estimated parameters associated with each 
%                                       model in data. numParams may be a scalar applied to all values in data, or a 
%                                       vector the same length as data. All elements of numParams must be 
%                                       positive integers.
%       numObs          -   [integer]   vector, a scalar or a numModels by 1 vector, sample sizes of the observed data 
%                                       associated with each models in data. numObs may be a scalar applied to all values in 
%                                       data, or a vector the same length as data. All elements numObs must be positive 
%                                       integers. 
%
%       where
%           numModels is the number of different to compare
%           numObs is the sample size
%           numVars is the number of variables include in multivariate models
%
% Optional Input Argument:
%
%       none
%
% Output Arguments:
%
%       results, a structure with the following fields:
%       infocrit.AIC	-   [double]	Akaike, numModels by 1 vector of statistics associated with each models
%       infocrit.AICc   -   [double]	Corrected Akaike, numModels by 1 vector of statistics associated with each models
%       infocrit.ABIC   -   [double]	Akaike-Bayesian, numModels by 1 vector of AIC statistics associated with each models
%       infocrit.SBIC	-   [double]	Schwarz-Bayesian, numModels by 1 vector of AIC statistics associated with each models
%       infocrit.HQIC	-   [double]	Hannan-Quinn, numModels by 1 vector of AIC statistics associated with each models
%       infocrit.FPE    -   [double]	Akaike Final Prediction Error, numModels by 1 vector of AIC statistics associated 
%                                       with each models
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       vec
%
% References:
%
%       Box, G.E.P., Jenkins, G.M., Reinsel, G.C., "Time Series Analysis: Forecasting and Control", 3rd edition, 
%       Prentice Hall, 1994.
%   
%       Liew, K. S. 2000. The performance of AICC as order determination criterion in the selection of ARMA time 
%       series models. Unpublished thesis, Department of Mathematics, Universiti Putra Malaysia.
%
%       Brockwell, P. J., Davis, R. A., 1996. Introduction to Time Series and Forecasting. Springer, New York.
%
%       Herman J. Bierens, Information Criteria and Model Selection, March 12, 2006 (lecture notes).
%
% Notes:
%       
%       TODO:   1)  Akaike Final Prediction Error (FPE) is not computed when mode = 'logL'. I need to find the 
%                   expression involving the log-likelihood.
%               2)  Akaike-Bayesian (ABIC) not yet implemented.
%               3)  Add more information criterion.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012
%           based on The MathWorks, Inc., 1999-2010 code "aicbic"


% Input and Output arguments checking
% ___________________________________

narginchk(4,4); nargoutchk(0,1);


% Compute information criterion
%______________________________

switch lower(mode)
	case 'logl'     % log-likelihood

        % Get the input more readable
        logL = data;
        
        % Input checking for mode = 'logL'
        % Ensure that data is a vector
        if numel(logL) == length(logL)      % check for a vector
            rowLL = (size(logL,1) == 1);	% flag a row vector for outputs
            if rowLL
                logL = vec(logL);           % convert to a column vector
            end
        else
            error('infocriterion:NonVectorlogL','infocriterion: data must be a vector when mode = ''logL'' ')
        end
        
        % Define the number of models
        numModels = size(logL,1);
        
        % Ensure numParams is a scalar, or compatible vector, of positive integers
        if numel(numParams) == length(numParams)	% check for a vector
            rowLL = (size(numParams,1) == 1);       % flag a row vector for outputs
            if rowLL
                numParams = vec(numParams);         % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumParams','infocriterion: numParams must be a vector')
        end

        if any(round(numParams) - numParams) || any(numParams <= 0)
            error('infocriterion:NonPositiveNumParams','infocriterion: numParams must be a positive integer')
        end

        if length(numParams) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numParams and data must be the same length.')
        end
        
        % Ensure numObs is a scalar, or compatible vector, of positive integers
        if numel(numObs) == length(numObs)	% check for a vector
            rowLL = (size(numObs,1) == 1);	% flag a row vector for outputs
            if rowLL
                numObs = vec(numObs);      % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumObs','infocriterion: NumObs must be a vector')
        end

        if any(round(numObs) - numObs) || any(numObs <= 0)
            error('infocriterion:NonPositiveNumObs','infocriterion: numObs must be a positive integer')
        end

        if length(numObs) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numObs and data must be the same length.')
        end
                
        % Compute information criteria
        % AIC
        infocrit.AIC = -2*logL + 2*numParams;
        % AICc
        infocrit.AICc = infocrit.AIC + 2*numParams.*(numParams+1)./(numObs-numParams-1);%2*((numParams+1).*(numParams+2))./(numObs-numParams-2);
        % ABIC - Not yet implemented
        infocrit.ABIC = [];
        % SBIC
        infocrit.SBIC = -2*logL + numParams.*log(numObs);
        % HQIC
        infocrit.HQIC = -2*logL + 2*numParams.*log(log(numObs));
        % FPE
        infocrit.FPE = [];
        
    case 'resuni'	% residuals from univariate models

        % Get the input more readable
        resid = data;
        
        % Input checking for mode = 'ResUni'
        % Define the number of models
        numModels = size(resid,1);
        
        % Ensure numParams is a scalar, or compatible vector, of positive integers
        if numel(numParams) == length(numParams)	% check for a vector
            rowLL = (size(numParams,1) == 1);       % flag a row vector for outputs
            if rowLL
                numParams = vec(numParams);         % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumParams','infocriterion: numParams must be a vector')
        end

        if any(round(numParams) - numParams) || any(numParams <= 0)
            error('infocriterion:NonPositiveNumParams','infocriterion: numParams must be a positive integer')
        end

        if length(numParams) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numParams and data must be the same length.')
        end
        
        % Ensure numObs is a scalar, or compatible vector, of positive integers
        if numel(numObs) == length(numObs)	% check for a vector
            rowLL = (size(numObs,1) == 1);	% flag a row vector for outputs
            if rowLL
                numObs = vec(numObs);       % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumObs','infocriterion: NumObs must be a vector')
        end

        if any(round(numObs) - numObs) || any(numObs <= 0)
            error('infocriterion:NonPositiveNumObs','infocriterion: numObs must be a positive integer')
        end

        if length(numObs) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numObs and data must be the same length.')
        end
        
        % Compute the estimated (MLE) variance
        sighat = zeros(numModels,1);
        for indexModel = 1:numModels 
            sighat(indexModel,:) = (resid(indexModel,:)*resid(indexModel,:)')/(numObs(indexModel));
        end
        
        % Compute information criteria
        % AIC
        infocrit.AIC  = log(sighat) + 2*numParams./numObs;
        % AICc
        infocrit.AICc = infocrit.AIC + 2*numParams.*(numParams+1)./(numObs-numParams-1);%2*((numParams+1).*(numParams+2))./(numObs-numParams-2);
        % SBIC
        infocrit.SBIC = log(sighat) + (numParams.*log(numObs))./numObs;
        % ABIC - Not yet implemented
        infocrit.ABIC = [];
        % HQIC
        infocrit.HQIC = log(sighat) + (2*numParams.*log(log(numObs)))./numObs;
        % FPE
        infocrit.FPE  = sighat.*(numObs+numParams)./(numObs-numParams);
        
    case 'resmulti'     % residuals from multivariate models

        % Get the input more readable
        resid = data;
        
        % Input checking for mode = 'ResMulti'
        % Define the number of models
        numModels = size(resid,1);
        
        % Define the number of endogenous variables
        numVars = size(resid,3);
        
        % Ensure numParams is a scalar, or compatible vector, of positive integers
        if numel(numParams) == length(numParams)	% check for a vector
            rowLL = (size(numParams,1) == 1);      % flag a row vector for outputs
            if rowLL
                numParams = vec(numParams);      % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumParams','infocriterion: numParams must be a vector')
        end

        if any(round(numParams) - numParams) || any(numParams <= 0)
            error('infocriterion:NonPositiveNumParams','infocriterion: numParams must be a positive integer')
        end

        if length(numParams) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numParams and data must be the same length.')
        end        

        % Ensure numObs is a scalar, or compatible vector, of positive integers
        if numel(numObs) == length(numObs)	% check for a vector
            rowLL = (size(numObs,1) == 1);	% flag a row vector for outputs
            if rowLL
                numObs = vec(numObs);      % convert to a column vector
            end
        else
            error('infocriterion:NonVectorNumObs','infocriterion: numObs must be a vector')
        end

        if any(round(numObs) - numObs) || any(numObs <= 0)
            error('infocriterion:NonPositiveNumObs','infocriterion: numObs must be a positive integer')
        end

        if length(numObs) ~= numModels
            error('infocriterion:VectorLengthMismatch', ...
                    'infocriterion: numObs and data must be the same length.')
        end
        
        % Compute the estimated (MLE) variance
        sighat = zeros(numModels,numVars,numVars);
        for indexModel = 1:numModels 
            sighat(indexModel,:,:) = (squeeze(resid(indexModel,:,:))'*squeeze(resid(indexModel,:,:)))/(numObs(indexModel));
        end
        
        % Compute information criteria
        % AIC
        infocrit.AIC = zeros(numModels,1);
        for indexModel = 1:numModels 
            infocrit.AIC(indexModel) = log(det(squeeze(sighat(indexModel,:,:)))) + 2*numParams(indexModel)/numObs(indexModel);
        end
        % AICc
        infocrit.AICc = infocrit.AIC + 2*numParams.*(numParams+1)./(numObs-numParams-1);%2*((numParams+1).*(numParams+2))./(numObs-numParams-2);
        % SBIC
        infocrit.SBIC = zeros(numModels,1);
        for indexModel = 1:numModels 
            infocrit.SBIC(indexModel) = log(det(squeeze(sighat(indexModel,:,:)))) + numParams(indexModel)*log(numObs(indexModel))/numObs(indexModel);
        end
        % ABIC - Not yet implemented
        infocrit.ABIC = [];
        % HQIC
        infocrit.HQIC = zeros(numModels,1);
        for indexModel = 1:numModels 
            infocrit.HQIC(indexModel) = log(det(squeeze(sighat(indexModel,:,:)))) + 2*numParams(indexModel)*log(log(numObs(indexModel)))/numObs(indexModel);
        end
        % FPE
        infocrit.FPE = zeros(numModels,1);
        for indexModel = 1:numModels 
            infocrit.FPE(indexModel) = det(squeeze(sighat(indexModel,:,:)))*((numObs(indexModel)+numParams(indexModel))/(numObs(indexModel)-numParams(indexModel)))^numVars;
        end
        
end % switch {mode}


end % function {infocriterion}

