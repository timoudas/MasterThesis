function [Y,checkRegimes,threshRatio] = tvarxnlsim(results,W,transVar,Y0,X)
% TVARXNLSIM: Simulates nonlinear (non-regime specific) sample paths of a multidimensional TVARX time series process.
%
% Syntax:
%
%       [Y,checkRegimes] = tvarxnlsim(results,W,transVar,Y0,X)
%
% Description:
%
%       TVARXNLSIM simulates nonlinear (non-regime specific) sample paths of a multidimensional TVARX time series process.
%
% Input Arguments:
%
%       results         -   [struct]            a results structure returned by TVARX function
%       W               -	[double]            numObs-by-numY-numPaths innovations data
%       transVar        -   [logical,double]    thVar is either:
%                                                   1 - a logical size(y,2)-by-1 vector indicating which variables (a single 
%                                                       value or a combination of endogenous and exogenous variables with same lag 
%                                                       order) to take for the transition variables.
%                                                           ex: if y is a 500-by-4, if thVar is logical([1 0 0 0]), 
%                                                           then the transition variable is the first one in y. 
%                                                   2 - a vector of exogenous (included in exo) or external (not included in y or exo) 
%                                                       transition variable (not included in y or exo). This vector has be already adjusted for the 
%                                                       time delay (results.timeDelay)
%
% Optional Input Arguments:
%
%       Y0              -	[double]            numObs-by-numPresampY0 or numObs-by-numPresampY0-by-numPaths presampled response data. Y0 is a  
%                                               matrix or a 3D array. If Y0 is empty or if numPresampY0 is less than the maximum AR lag in  
%                                               results, presample values are padded with zeros. If numPresampY0 is greater than the maximum AR  
%                                               lag, the most recent samples from the last rows of each path of Y0 are used. If Y has multiple 
%                                               paths, Y0 must contain either a single path (applied to all paths in Y) or at least as many 
%                                               paths as in Y (extra paths are ignored).
%                                                   Default: Y0 = []
%       X               -	[double]            numObs-by-numX or numObs-by-numX-by-numPaths exogenous data. X is a matrix or a 3D array. If 
%                                               Y has multiple paths, X must contain either a single path (applied to all paths in Y) or at 
%                                               least as many paths as in y (extra paths are ignored).
%                                                   Default: X = []
%
% Output Arguments:
%
%       Y               -	[double]            numObs-by-numVars-by-numPaths array of simulated paths
%       checkRegimes	-	[double]            numObs-by-numPaths array of regimes indices
%       threshRatio     -	[double]            numThresh-by-nnumPaths array of ratio of period passed into a specific regime
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions 
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

narginchk(4,5); nargoutchk(0,3);

[results,W,transVar,Y0] = checkInput(results,W,transVar,Y0);

% Define constant
numRegimes = size(results,1);
[numObs,~,numPaths] = size(W);
numLags = results(1).sizes.numLags;
numThresh = results(1).sizes.numThresh;


% Options
% _______

if nargin < 5 
	X = double.empty(size(W,1),0);
end

X = checkOptions(X,results,W);


% Cholesky decomposition of the residuals variance-covariance matrix
% __________________________________________________________________

sigechol = cell(numRegimes,1);
for indRegimes = 1:numRegimes
    sigechol{indRegimes} = chol(results(indRegimes).stats.sige);
end


% Initialization
% ______________

Y            = zeros(size(W));
checkRegimes = zeros(size(W,1),numPaths);
ARlag        = 1:numLags;
nAR          = max(ARlag);
timeDelay    = results(1).transition.timeDelay;
thresh       = results(1).threshold.thresh;


isConstant = results(1).constant;
isTrend    = results(1).trend;
isX        = (~isempty(results(1).beta.exo) && ~isempty(X));


% Simulate data
% _____________

for indPaths = 1:numPaths
	for indObs = 1:numObs

        % Determine the transition variable
        if islogical(transVar)
            if indObs <= timeDelay
            	transVarLevel = Y0(nAR-timeDelay+indObs,transVar,indPaths);
            else
            	transVarLevel = Y(indObs-timeDelay,transVar,indPaths);
            end
        elseif isnumeric(transVar)
        	transVarLevel = transVar(indObs);
        end
        
        % Check the threshold to determine in which regime we are
        switch numThresh
        	case 1  % One threshold
                if transVarLevel <= thresh
                    indRegimes = 1;
                elseif transVarLevel > thresh
                	indRegimes = 2;
                end
            case 2  % Two thresholds
                if transVarLevel <= thresh(1)
                	indRegimes = 1;
                elseif transVarLevel > thresh(1) && transVarLevel <= thresh(2)
                	indRegimes = 2;
                elseif transVarLevel > thresh(2)
                	indRegimes = 3;
                end
        end
        checkRegimes(indObs,indPaths) = indRegimes;
        
        % Add constant term to response:
        if isConstant
			Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + results(indRegimes).beta.cons';
        end

        % Add linear trend term to response:
        if isTrend
			Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + indObs*results(indRegimes).beta.trend';
        end
        
		% Add exogenous inputs to response:
        if isX
            Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (X(indObs,:,indPaths)*results(indRegimes).beta.exo');
        end        

        % Add autoregressive terms to response:
        for indLags = 1:numLags
            if indObs <= ARlag(indLags)
				Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (Y0(nAR-ARlag(indLags)+indObs,:,indPaths)*results(indRegimes).beta.lags{indLags}');
            else
                Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (Y(indObs-ARlag(indLags),:,indPaths)*results(indRegimes).beta.lags{indLags}');
            end
        end % for {indLags}
        
        % Add current innovation to response:
        Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + W(indObs,:,indPaths)*sigechol{indRegimes};

	end % for {indObs}
end % for {indPaths}

% Compute the threshold period over total period
threshRatio = NaN(numRegimes,numPaths);
for indPaths = 1:numPaths
    for indRegimes = 1:numRegimes
        threshRatio(indRegimes,indPaths) = (sum(double(checkRegimes(checkRegimes(:,indPaths) == indRegimes,indPaths)))/indRegimes)/length(checkRegimes(:,indPaths));
    end
end


end % function {tvarxnlsim}


% ---------------------------------
function [results,W,transVar,Y0] = checkInput(results,W,transVar,Y0)
% checkInput: Local function to check the validity of required inputs

% Set constant
numRegimes = size(results,1);
[~,numDims,numPaths] = size(W); % Define size
numLags = results(1).sizes.numLags;

% results
validateattributes(results,{'struct'},{},'tvarxnlsim','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxnlsim:InputError','tvarxnlsim: The results structure provided must be returned by TVARX function')
end

% W
for indRegimes = 1:numRegimes
    validateattributes(W,{'numeric'},{'real','finite','size',[NaN numDims NaN]},'tvarxnlsim','W',2);
end

% transVar
validateattributes(transVar,{'logical','numeric'},{},'tvarxnlsim','transVar',3);
if islogical(transVar)                                                                                              
    validateattributes(transVar,{'logical'},{'numel',size(W,2)},'tvarxnlsim','transVar',3);
    if sum(double(transVar)) > 1
        error('tvarxnlsim:InputError','tvarxnlsim: currently, for transVar, only a single value is implemented')
    elseif sum(double(transVar)) == 0
        error('tvarxnlsim:InputError','tvarxnlsim: at least one transition variable must be declared')
    end
elseif isnumeric(transVar)
    validateattributes(transVar,{'numeric'},{'real','finite','nrows',size(W,1),'ncols',1},'tvarxnlsim','transVar',3);
end

% Y0
validateattributes(Y0,{'numeric'},{'real','finite','nrows',numLags,'ncols',numDims},'tvarxnlsim','Y0',4);
if size(Y0,3) ~= 1 && size(Y0,3) ~= numPaths
    error('tvarxnlsim:InputError','tvarxnlsim: the third dimension size of Y0 should be 1 or equal to its W counterpart.');
end
if size(Y0,3) == 1
	Y0 = repmat(Y0,[1 1 numPaths]);
end


end % subfunction {checkInput}


% ---------------------------------
function X = checkOptions(X,results,W)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Set constant
[numObs,~,numPaths] = size(W); % Define size

% X
if ~isempty(X)
    validateattributes(X,{'numeric'},{'real','finite','nrows',numObs,'ncols',size(results(1).beta.exo,2)},'tvarxnlsim','X',5);
    if size(X,3) ~= 1 && size(X,3) ~= numPaths
        error('tvarxnlsim:InputError','tvarxnlsim: the third dimension size of X should be 1 or equal to its W counterpart.');
    end
    if size(X,3) == 1
        X = repmat(X,[1 1 numPaths]);
    end
end


end % subfunction {checkOptions}

