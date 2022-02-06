function results = tvarx(y,numLags,numThresh,transVar,timeDelay,options)
% TVARX: Estimate an unrestricted TVARX model by ordinary least squares or a restricted TVARX model by generalized least squares.
%
% Syntax:
%
%       results = tvarx(y,numLags,numThresh,transVar,timeDelay)
%       results = tvarx(y,numLags,numThresh,transVar,timeDelay,options)
%       results = tvarx(y,numLags,numThresh,transVar,timeDelay,struct('cons',true,'trend',true))
%
% Description:
%
%       Estimate a threshold vector autoregressive with (potentially) exogenous variables model by least squares of the following form:
%
%       y(t) = [A11*y(t-1) + ... + A1p*y(t-p) + B1*x(t) + C1*D(t) + e1(t)] * 1(s(t-d) <= gamma) + 
%              [A21*y(t-1) + ... + A2p*y(t-p) + B2*x(t) + C2*D(t) + e2(t)] * 1(s(t-d) > gamma)
%
%       The unrestricted TVARX model is estimated by sequential conditional ordinary least squares (equivalent to maximum likelihood) and the restricted TVARX 
%       model is estimated by sequential conditional generalized least squares (equivalent to restricted maximum likelihood).
%
%       TVARX treats NaNs in 'y' and 'exo' as missing values, and removes them.
%
%       TVARX can automatically include a constant and a linear trend. Other deterministic regressors must be included into the matrix of exogenous regressors 
%       'exo'.
%
% Input Arguments:
%
%       y               -	[double]            numObs-by-numY, dependent variable matrix (or response observations, or regressands)
%       numLags         -	[integer]           scalar, lag length (common to all regimes)
%       numThresh       -	[integer]           scalar, number of thresholds. Currently, the number of thresholds could be 1 (= 2 regimes)
%                                               or 2 (= 3 regimes).
%       transVar        -	[logical,double]	transVar is either:
%                                                   1 - a logical numY-by-1 vector indicating which variables (a single value or a combination of endogenous 
%                                                       and exogenous variables with same lag order) to take for the transition variables.
%                                                       ex: if y is a 500-by-4 , if transVar is logical([1 0 0 0]), then the transition variable is the first 
%                                                           one in y. 
%                                                   2 - a numObs-by-1 vector of exogenous (included in exo) or external (not included 
%                                                       in y or exo) transition variable.
%       timeDelay       -   [integer]           non-negative integer, 'time delay' for the threshold variable. 
%                                                   If transVar is logical, an endogenous variables is the transition variables and timeDelay is 
%                                                       not lower than 1 and not greater than numLags. 
%                                                   If transVar is numeric, an exogenous or external variables is the transition variables and 
%                                                       timeDelay is not lower than 0 and not greater than numLags. 
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values (applied to all regimes, separately):
%       Model Specification options:
%       	cons        -   [logical]	a constant
%                                           true	include a constant  
%                                           false	does not include a constant (default)
%       	trend       -   [logical]	a linear trend
%                                           true	include a linear trend
%                                           false	does not include a linear trend (default)
%       	exo         -	[double]	numObs-by-numVars, independent variables matrix (regressors). Its a design matrix, with rows corresponding to 
%                                       observations and columns to predictor variables graph
%                                       	Default: empty
%       	ylog        -   [logical]   a logical numY-by-1 vector, log the endogenous series before estimating 
%                                       the VARX: true means apply the log, false means no log is taken 
%       	                            	Default: No log (logical(zeros(numY,1)))	
%       	ydiff       -   [logical]   a logical numY-by-1 vector, first difference the endogenous series before estimating 
%                                       the VARX: true means apply the difference, false means no difference is taken 
%       	                            	Default: No difference (logical(zeros(numY,1)))
%           Rest        -	[struct]	structure of numerical restrictions. All restrictions are defined this way:
%                                           1 - NaN means no restrictions applied on this parameters, or
%                                           2 - a real numeric value to apply this numerical restriction on a specific parameter
%                                       The structure has these specifics fields:
%                                       	ARsolve     a numY-by-numY-by-numLags matrix (default: NaN(numY,numY,numLags))
%                                                       if more than numLags layer are provided, only the first numLags will be kept
%                                           Csolve      a numY-by-1 matrix (default: if cons = true, a NaN(numY,1))
%                                           Tsolve      a numY-by-1 matrix (default: if trend = true, a NaN(numY,1))
%                                           EXOsolve	a numY-by-numX matrix (default: if exo not empty, a NaN(numY,numX))
%           algo        -   [char]      string, way to compute the MLE estimator of beta:
%                                           'explicit'  explicit formula for global minimum
%                                           'qr'        QR decomposition (default)
%                                           'svd'       SVD decomposition
%       Variance-Covariance options:
%       	commonvar	-   [logical]	each regime share the same variance-covariance matrix
%                                           true    same distribution
%                                           false	specific distribution (default)
%       	het         -   [logical]	type of covariance estimator
%                                           true    heteroskedastic
%                                           false	homoskedastic (default)
%       	corre       -   [logical]	the assumed structure of the error covariance matrix
%                                           true    correlated residuals
%                                           false	uncorrelated residuals (default)
%       Threshold options:
%       	optimcrit   -   [char]      string, criterion to optimize
%                                           'loglik'    for log-likelihood (default)
%                                           'logdet'    for log of the determinant of the residuals variance-covariance matrix
%                                           'ssr'       for sum of square residuals
%       	numGrid     -   [integer]   number of equally spaced elements of the grid (especially useful for numThresh = 2)
%                                           Default: 50 (for each threshold)
%       	gamm        -	[double]    vector, prespecified threshold values. No estimation is done when gamm is not empty
%                                           ex: options.gamm = 50 (if numThresh = 1)
%                                               options.gamm = [25 50] (if numThresh = 2)
%                                           Default: empty
%       	trim        -   [double]    scalar, trimming parameter indicating the minimal percentage of observations in each 
%                                       regime.
%                                       	Default: 15 (for 15%)
%           around      -	[logical]	perform a second optimization round around the optimal threshold found in the first
%                                       round (Note: You can also use more points in the grid instead of this option.)
%                                       	true    perform a second optimization round
%                                           false	do not perform a second optimization round (default)
%       Results Presentation options:
%       	prtr        -	[logical]	print results
%       	                                true	print results
%       	                                false	do not print results (default)
%       	vnames      -   [char]      string, vector of variable names (only relevant when prtr = true). If 
%       	                            prtr = true and no vnames is provided, a string when default variables names will
%                                       generated
%       	                            	e.g.	vnames = char('y1','y2','x1','x2') 
%       	graph       -	[logical]	graph fitted vs observed values
%       	                                true	graph fitted vs observed values
%       	                                false	do not graph fitted vs observed values (default)
%           calstruct	-   [struct]	a calendar struct, as returned by cal.m
%
% Output Arguments:
%
%       results, a structure with the following fields:         
%       Common to all threshold:
%           meth        -   [char]      method
%           description	-   [char]      method description
%           fitmeth     -   [char]      fitting method
%           algo        -   [char]      optimization algorithm
%           optim       -   [struct]    structure, whose fields are optimization information
%                                       	optimcrit       -   [char]      string, criterion to optimize
%                                           gridpt          -   [double]    numGrid-by-numThresh, matrix of potential threshold
%                                       If optimcrit = 'loglik':
%                                       	logL            -   [double]	scalar, log-likelihood evaluated at beta
%                                           logLloop        -   [double]	numGrid-by-numThresh, matrix of log-likelihood for each 
%                                                                           threshold
%                                       If optimcrit = 'logdet':
%                                       	logDet          -   [double]	scalar, log of the determinant of the residuals variance
%                                                                           -covariance matrix evaluated at beta
%                                           logDetloop      -   [double]	numGrid-by-numThresh, matrix of log of the determinant 
%                                                                           of the residuals variance-covariance matrix
%                                       If optimcrit = 'ssr':
%                                       	ssr             -   [double]	scalar, sum of square residuals evaluated at beta
%                                           ssrloop         -   [double]	numGrid-by-numThresh, matrix of sum of square residuals
%           threshold	-   [struct]    structure, whose fields are threshold information
%                                           thresh          -	[double]	numThresh-by-1, thresholds value(s)
%                                           gamm            -   [double]    vector, prespecified threshold values. No estimation is done 
%                                                                           when gamm is not empty 
%           transition	-   [struct]    structure, whose fields are transition variable information
%                                       	transVar        -   [logical,double]	logical or double, the input to the function
%                                           transVarLevel	-   [double]	numObs-by-1, transition variables
%                                           transVarInd     -   [double]	numObs-by-1, vector for the position of regimes
%                                                                               If numRegimes = 2, take the value 1 and 2    
%                                                                               If numRegimes = 3, take the value 1, 2 and 3   
%                                           timeDelay       -   [double]    scalar, 'time delay' for the threshold variable 
%                                           trim            -   [integer]	scalar, trimming parameter indicating the minimal percentage of 
%                                                                           observations in each regime 
%           sizes       -   [struct]    structure, whose fields are defined size
%                                           numGrid         -   [double]	scalar, number of points in the grid 
%                                           numLags         -   [double]    scalar, lag length
%                                           numThresh       -   [double]    scalar, number of threshold
%                                           numRegimes      -   [double]    scalar, number of regime (= numThresh+1)
%           constant    -   [logical]	logical, true if constant
%           trend       -   [logical]   logical, true if linear trend
%           ylog        -   [logical]	logical numY-by-1 vector, log the endogenous series before estimating 
%           ydiff       -   [logical]	logical numY-by-1 vector, first difference the endogenous series before estimating 
%           commonvar	-   [logical]   logical, true if all regimes share the same variance-covariance matrix  
%           het         -   [logical]   logical, true if residuals are assumed to be hetetoskedastic   
%           corre       -   [logical]   logical, true if residuals are assumed to be correlated
%           prtr        -   [logical]   true if print results
%           graph       -   [logical]   true if graph fitted values vs observed values
%           vnames      -   [char]      string, vector of variable names (only relevant when prtr = true).
%           calstruct   -   [struct]    structure, a calendar structure
%           options     -   [struct]	structure, options provided by the user
%       Specific to each threshold:
%           solmat      -   [struct]	structure, whose fields are matrices of the solution
%                                       If algo = 'explicit':
%                                       	nothing
%                                       If algo = 'qr':
%                                       	Q               - 	[double]    (numObs-numNaN)-by-numX, Q from the QR decomposition
%                                                                           of the design matrix
%                                           R               - 	[double]    numX-by-numX, R from the QR Decomposition of the
%                                                                   design matrix    
%                                           perm            -   [double]    1-by-numX, permutation information
%                                       If algo = 'svd':
%                                       	U               - 	[double]    (numObs-numNaN)-by-(numObs-numNaN), U from the SVD 
%                                                                           decomposition of the design matrix    
%                                           V               - 	[double]    numX-by-numX, V from the SVD decomposition of the 
%                                                                           design matrix
%                                           S               -   [double]    (numObs-numNaN)-by-numX, rectangular diagonal matrix
%                                                                           with nonnegative real numbers on the diagonal
%           logCL       -   [double]	numObs-by-1, vector of log-likelihoods corresponding to beta
%           threshold	-   [struct]    structure, whose fields are threshold information
%                                           threshRatio     -   [double]	scalar, proportion of observations in a specific regime
%           data        -	[struct]	structure, whose fields are matrices are the data
%                                           y               -   [double]    numObs-by-numY, endogenous variables
%                                           exo             -   [double]    numObs-by-numExo, exogenous variables
%                                           yreg            -   [double]    (numObs-numLags)-by-numY, regressands (endogenous)
%                                           xreg            -   [double]    (numObs-numLags)-by-numX, regressors (both endogenous and exogenous)
%                                           ynan            -   [double]    (numObs-numNaN)-by-numY, regressands with NaN removed
%                                           xnan            -   [double]	(numObs-numNaN)-by-numX, regressors with NaN removed
%           beta        -   [struct]	structure, estimated coefficients
%           stderr      -   [struct]    structure, standard errors of estimated coefficients
%           tstat       -   [struct]    structure, t-stats of estimated coefficients
%           pval        -   [struct]    structure, p-values of estimated coefficients
%           rest        -   [struct]    structure, restrictions of estimated coefficients
%           yhat        -   [double]	(numObs-numNaN)-by-numY, fitted values
%           resid       -   [double]    (numObs-numNaN)-by-numY, residuals
%           stats       -   [struct]    structure, whose fields are matrices are the statistics
%                                       	r2              -   [double]    numY-by-1, R-squared
%                                           rbar2           -   [double]    numY-by-1, adjusted R-squared
%                                           vcov            -   [double]    numX-by-numX, parameters' variance-covariance matrix 
%                                           sigu            -   [double]    numY-by-numY, sum of squared residuals    
%                                           sige            -   [double]    numY-by-numY, residuals variance    
%           sizes       -   [struct]    structure, whose fields are defined size
%                                           numObs          -   [double]	scalar, number of observations
%                                           numObsLag       -   [double]    scalar, number of observations adjusted to feed the lags
%                                           numObsNaNLag    -   [double]    scalar, number of observations adjusted to feed the lags after NaN are removed
%                                           numY            -   [double]    scalar, number of regressands
%                                           numX            -   [double]    scalar, total number of regressors
%                                           numExo          -   [double]    scalar, number of exogenous regressors
%                                           numNaN          -   [double]    scalar, number of removed values (NaN)
%                                           numParams       -   [double]    scalar, number of parameters
%                                           numActive       -   [double]    scalar, number of active (unrestricted) parameters (including variance-covariance matrix)
%           allNaN      -   [logical]	logical, true if the series is all NaN
%           rankdef     -   [logical]   logical, true if the series is rank deficient    
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, lagmat, mvlinreg, demean, mdmax, mdmin, tvarxloglik, varxvcov, invprctile, tvarxprt, tvarxplt
%
% References:
%
%       none
%
% Notes:
%
%       TODO:	1) Combination of variables for the transition variable (currently only one variable can be used as threshold).
%               2) Include DiRect, Pattern Search, Simulated Annealing and Genetic algorithms as optimization option (currently, only the 
%                  Grid Search algorithm is implemented).
%               3) Same constant and trend. Currently, only regime specific constant and trend can be estimated.
%               4) Regime specific restrictions. Currently, only the same set of restrictions can be estimated.
%
%       Differences:	1)	A difference in the estimation of the parameters with other VARX package when lag is > 2 is possible. 
%                           This is because I initialize the trend with the number of lags (i.e., when numLags = 2, the trend is [2 3 ...numObs]), 
%                           while other package may initialize the trend with 1. 
%
%       Known bugs:	1)	Do not use ylog and ydiff for the moment. The other functions are not yet updated to reflect these new options.
%                   2)  The options "around" is not working very well. Do not use this option for the moment. Increase the number of grid points instead.
%                   3)  Alternative estimation specifications to be possibbly included in the future:
%                       The first one estimate this model:
%                           y(t) = [A11*y(t-1) + ... + A1p*y(t-p) + B1*x(t) + C1*D(t) + e1(t)]
%                                  [A21*y(t-1) + ... + A2p*y(t-p) + B2*x(t) + C2*D(t) + e2(t)] * 1(s(t-d) > gamma)
%                       with this code:
%                           res = mvlinreg({yreg},...
%                                   {[xreg xreg.*repmat(double(transVarLevel > threshGrid{1}(iter)),1,size(xreg,2))]},'ols',lower(algo));
%                       The second one allows to estimate common constant, trend, and exo variables:
%                           y(t) = B*x(t) + C*D(t) +  
%                                  [A11*y(t-1) + ... + A1p*y(t-p) + e1(t)] * 1(s(t-d) <= gamma) + 
%                                  [A21*y(t-1) + ... + A2p*y(t-p) + e2(t)] * 1(s(t-d) > gamma)
%                       with this code:
%                           res = mvlinreg({yreg},...
%                                   {[xreg.*repmat(double(transVarLevel <= threshGrid{1}(iter)),1,size(xreg,2)) ...
%                                   xreg.*repmat(double(transVarLevel > threshGrid{1}(iter)),1,size(xreg,2))]},'ols',lower(algo));
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(5,6); nargoutchk(1,1);

[y,numLags,numThresh,transVar,timeDelay] = checkInput(y,numLags,numThresh,transVar,timeDelay);


% Options
% _______

if nargin < 6
    options = struct.empty;
elseif nargin == 6
	validateattributes(options,{'struct'},{},'tvarx','options',5);
end

[cons,trend,exo,ylog,ydiff,Rest,algo,commonvar,het,corre,optimcrit,numGrid,gamm,trim,around,prtr,graph,vnames,calstruct] = checkOptions(options,y,numThresh);


% Restrictions
% ____________

Restrictions = double.empty(size(y,2),0);  % Initialize restrictions matrix

% Get default or user-provided restrictions
getOptions(Rest, ...
	'ARsolve',	NaN(size(y,2),size(y,2),numLags), ...
    'Csolve',	NaN(size(y,2),double(cons)), ...
    'Tsolve',	NaN(size(y,2),double(trend)), ...
    'EXOsolve',	NaN(size(y,2),size(exo,2)));

% AR solve
if size(ARsolve,3) <= numLags  %#ok<NODEF>	
    validateattributes(ARsolve,{'numeric'},{'size',[size(y,2) size(y,2) numLags]},'tvarx','options.Rest.ARsolve',6);
elseif size(ARsolve,3) > numLags
    ARsolve = ARsolve(:,:,1:numLags);
    validateattributes(ARsolve,{'numeric'},{'size',[size(y,2) size(y,2) numLags]},'tvarx','options.Rest.ARsolve',6);
end

% Csolve
validateattributes(Csolve,{'numeric'},{'size',[size(y,2) double(cons)]},'tvarx','options.Rest.Csolve',6);

% Tsolve
validateattributes(Tsolve,{'numeric'},{'size',[size(y,2) double(trend)]},'tvarx','options.Rest.Tsolve',6);

% EXOsolve
validateattributes(EXOsolve,{'numeric'},{'size',[size(y,2) size(exo,2)]},'tvarx','options.Rest.EXOsolve',6);


% Check inputs, form inputs and restrictions matrices 
% ___________________________________________________

% Check endogenous inputs
numObs = size(y,1); % Define numObs, the number of observations
Restrictions = [Restrictions reshape(ARsolve,[size(y,2) size(y,2)*numLags 1])];

% Check exogenous inputs
switch ~isempty(exo)
    case 0
        % No exogenous regressors
    case 1
        Restrictions = [Restrictions EXOsolve];
end

% Check constant
switch cons
    case 0
        % No constant
    case 1
        exo = [exo ones(numObs,1)];
        Restrictions = [Restrictions Csolve];
end

% Check trend
switch trend
    case 0
        % No linear trend
    case 1
        exo = [exo (1:1:numObs)'];
        Restrictions = [Restrictions Tsolve];
end

% Check final exogenous inputs
numExo = size(exo,2);

% Form the regressands and regressors matrices
switch ~isempty(exo)
	case 0  % No exogenous variables
        yreg = y(numLags+1:numObs,:);
        xreg = lagmat(y,1:1:numLags);
        xreg = xreg(numLags+1:numObs,:);
	case 1  % With exogenous variables (constant, linear trend, and/or others)
        yreg = y(numLags+1:numObs,:);
        xreg = [lagmat(y,1:1:numLags) exo];
        xreg = xreg(numLags+1:numObs,:);
end   

% Form the final restrictions matrices, if any
isRest = all(all(isnan(Restrictions)));
switch isRest
	case 1	% No restrictions
        % Nothing
    case 0  % Restrictions
		Rest = eye(size(yreg,2)*size(xreg,2));
		r = Restrictions(:);
	
		Rest = Rest(any(isnan(r),2),:).';
		r(isnan(r)) = 0;
end


% Form the threshold vector
% _________________________

if islogical(transVar)
    transVarLevel = y(:,transVar);
elseif isnumeric(transVar)
    transVarLevel = transVar;
end
transVarLevel = lagmat(transVarLevel,timeDelay);
transVarLevel = transVarLevel(numLags+1:numObs);
threshSort    = sort(transVarLevel);


% Optimization loop: Find the threshold 
% _____________________________________

switch ~isempty(gamm)
	case 1	% Threshold provided by user

        % Initialize array (for structure results only)
        logLloop = NaN;
        logDetloop = NaN;

        % Identify threshold value(s)
        gridpt = gamm;
        thresh = gamm;

    case 0	% Threshold not provided by user, need to be estimated

        % Check "around" for a second optimization round
        for i = 1:double(around)+1
            
            % Build grid point
            switch i
            	case 1
                    gridpt = linspace(prctile(transVarLevel,trim),prctile(transVarLevel,100-trim),numGrid)';
                    threshGrid = cell(numThresh,1);
                    [threshGrid{:}] = ndgrid(gridpt);
                case 2
                    threshGrid = cell(numThresh,1);
                    p = invprctile(threshSort,thresh);
                    checkPrctile = (p - 1);
                    switch numThresh
                        case 1  % One threshold
                            gridpt = linspace(prctile(transVarLevel,checkPrctile-trim/2),prctile(transVarLevel,checkPrctile+trim/2),numGrid)';
                        case 2  % Two thresholds
                            gridpt = linspace(prctile(transVarLevel,checkPrctile(1)-trim/2),prctile(transVarLevel,checkPrctile(2)+trim/2),numGrid)';
                    end
                    [threshGrid{:}] = ndgrid(gridpt);
            end % switch {i}
            
            % Initialize array
            logLloop   = NaN(size(threshGrid{1}));
            logDetloop = NaN(size(threshGrid{1}));
            ssrloop    = NaN(size(threshGrid{1}));

            switch numThresh
                case 1  % One threshold

                    % Find the threshold: Optimization loop
                    for iter = 1:numGrid

                        % Estimates both low and high regimes
                        switch isRest
                            case 1	% No restrictions
                                res = mvlinreg({yreg(transVarLevel <= threshGrid{1}(iter),:) ...
                                                yreg(transVarLevel >  threshGrid{1}(iter),:)}',...
                                               {xreg(transVarLevel <= threshGrid{1}(iter),:) ...
                                                xreg(transVarLevel >  threshGrid{1}(iter),:)}','ols',lower(algo));
                            case 0	% Restrictions
                                res = mvlinreg({yreg(transVarLevel <= threshGrid{1}(iter),:) ...
                                                yreg(transVarLevel >  threshGrid{1}(iter),:)}',...
                                               {xreg(transVarLevel <= threshGrid{1}(iter),:) ...
                                                xreg(transVarLevel >  threshGrid{1}(iter),:)}','rls',lower(algo), ...
                                               {Rest Rest}',{r r}');
                        end % switch {checkRest}

                        % Compute the residuals
                        resid = cell(size(res,1),1);
                        for indRegimes = 1:size(res,1)
                            resid{indRegimes} = res(indRegimes).data.ynan - res(indRegimes).data.xnan*res(indRegimes).beta;
                        end

                        % Optimization criterion
                        switch lower(optimcrit)
                            case 'loglik'	% Compute the log-likelihood
                                logLloop(iter) = tvarxloglik(resid,commonvar);

                            case 'logdet'   % Compute the log of the determinant of the residuals variance-covariance matrix
                                logDetloop(iter) = tvarxlogdet(resid,commonvar);

                            case 'ssr'      % Compute the sum of square residuals
                                residall = [];
                                for indRegimes = 1:size(resid,1)
                                    residall = [residall;resid{indRegimes}]; %#ok<AGROW>
                                end
                                ssrloop(iter) = sum(residall(:).^2);

                        end % switch {optimcrit}
                    end % for {iter}

                case 2  % Two thresholds
                    
                    % Find the threshold: Optimization loop
                    for iter1 = 1:numGrid
                        for iter2 = 1:numGrid

                            % Check if the grid points respect the trim parameter
                            p = invprctile(threshSort,[threshGrid{1}(iter1,iter2) threshGrid{2}(iter1,iter2)]);
                            checkTrim = ((p(2) - p(1)) >= trim);
                            
                            switch checkTrim
                                case 0  % Does not respect the trim parameter
                                    switch lower(optimcrit)
                                        case 'loglik'	% Maximize the log-likelihood
                                            logLloop(iter1,iter2) = -Inf;
                                        
                                        case 'logdet'   % Maximize the log of the determinant of the residuals variance-covariance matrix
                                            logDetloop(iter1,iter2) = Inf;
                                    
                                        case 'ssr'      % Maximize the sum of square residuals
                                            ssrloop(iter1,iter2) = Inf;
                                    
                                    end % switch {optimcrit}

                                case 1  % Respect the trim parameter
                                    % Estimates low, middle and high regimes
                                    switch isRest
                                        case 1	% No restrictions
                                            res = mvlinreg({yreg(transVarLevel <= threshGrid{1}(iter1,iter2),:) ...
                                                            yreg(transVarLevel >  threshGrid{1}(iter1,iter2) & transVarLevel <= threshGrid{2}(iter1,iter2),:) ...
                                                            yreg(transVarLevel >  threshGrid{2}(iter1,iter2),:)}',...
                                                           {xreg(transVarLevel <= threshGrid{1}(iter1,iter2),:) ...
                                                            xreg(transVarLevel >  threshGrid{1}(iter1,iter2) & transVarLevel <= threshGrid{2}(iter1,iter2),:) ...
                                                            xreg(transVarLevel >  threshGrid{2}(iter1,iter2),:)}','ols',lower(algo));
                                        case 0	% Restrictions
                                            res = mvlinreg({yreg(transVarLevel <= threshGrid{1}(iter1,iter2),:) ...
                                                            yreg(transVarLevel >  threshGrid{1}(iter1,iter2) & transVarLevel <= threshGrid{2}(iter1,iter2),:) ...
                                                            yreg(transVarLevel >  threshGrid{2}(iter1,iter2),:)}',...
                                                           {xreg(transVarLevel <= threshGrid{1}(iter1,iter2),:) ...
                                                            xreg(transVarLevel >  threshGrid{1}(iter1,iter2) & transVarLevel <= threshGrid{2}(iter1,iter2),:) ...
                                                            xreg(transVarLevel >  threshGrid{2}(iter1,iter2),:)}','rls',lower(algo),...
                                                            {Rest Rest Rest}',{r r r}');
                                    end % switch {checkRest}

                                    % Optimization criterion
                                    resid = cell(size(res,1),1);
                                    for indRegimes = 1:size(res,1)
                                        resid{indRegimes} = res(indRegimes).data.ynan - res(indRegimes).data.xnan*res(indRegimes).beta;
                                    end

                                    switch lower(optimcrit)
                                        case 'loglik'	% Compute the log-likelihood
                                            logLloop(iter1,iter2) = tvarxloglik(resid,commonvar);

                                        case 'logdet'   % Compute the log of the determinant of the residuals variance-covariance matrix
                                            logDetloop(iter1,iter2) = tvarxlogdet(resid,commonvar);

                                        case 'ssr'      % Compute the sum of square residuals
                                            residall = [];
                                            for indRegimes = 1:size(resid,1)
                                                residall = [residall;resid{indRegimes}]; %#ok<AGROW>
                                            end
                                            ssrloop(iter1,iter2) = sum(residall(:).^2);

                                    end % switch {optimcrit}

                            end % switch {checkTrim}
                        end % for {iter2}
                    end % for {iter1}

            end % switch {numThresh}
                    
            % Find the position of the threshold that maximized the log-likelihood
            switch lower(optimcrit)
                case 'loglik'	% Maximize the log-likelihood
                	[~,Icrit] = mdmax(logLloop);

                case 'logdet'   % Minimize the log of the determinant of the residuals variance-covariance matrix
                    [~,Icrit] = mdmin(logDetloop);
            
                case 'ssr'      % Minimize % Compute the sum of square residuals
                    [~,Icrit] = mdmin(ssrloop);

            end % switch {optimcrit}

            % Identify threshold value
            switch numThresh
                case 1  % One threshold                
                    thresh = threshGrid{1}(Icrit(1),Icrit(2));
                
                case 2  % Two thresholds                    
                    thresh = [threshGrid{1}(Icrit(1),Icrit(2)) threshGrid{2}(Icrit(1),Icrit(2))];
            
            end % switch {numThresh}
            
        end % for {i}
end % switch {checkGamma}


% Estimation
% __________

switch numThresh
    case 1  % One threshold
        % Estimates both low and high regimes
        switch isRest
            case 1	% No restrictions
            	res = mvlinreg({yreg(transVarLevel <= thresh,:) yreg(transVarLevel > thresh,:)}',...
                               {xreg(transVarLevel <= thresh,:) xreg(transVarLevel > thresh,:)}','ols',lower(algo));
            
            case 0	% Restrictions
            	res = mvlinreg({yreg(transVarLevel <= thresh,:) yreg(transVarLevel > thresh,:)}',...
                               {xreg(transVarLevel <= thresh,:) xreg(transVarLevel > thresh,:)}','rls',lower(algo),...
                               {Rest Rest}',{r r}');
        
        end % switch {checkRest}
            
    case 2  % Two thresholds
        % Estimates both low and high regimes
        switch isRest
            case 1	% No restrictions
                 res = mvlinreg({yreg(transVarLevel <= thresh(1),:) ...
                                 yreg(transVarLevel >  thresh(1) & transVarLevel <= thresh(2),:) ...
                                 yreg(transVarLevel >  thresh(2),:)}',...
                                {xreg(transVarLevel <= thresh(1),:) ...
                                 xreg(transVarLevel >  thresh(1) & transVarLevel <= thresh(2),:) ...
                                 xreg(transVarLevel >  thresh(2),:)}','ols',lower(algo));
             
            case 0	% Restrictions
                 res = mvlinreg({yreg(transVarLevel <= thresh(1),:) ...
                                 yreg(transVarLevel >  thresh(1) & transVarLevel <= thresh(2),:) ...
                                 yreg(transVarLevel >  thresh(2),:)}',...
                                {xreg(transVarLevel <= thresh(1),:) ...
                                 xreg(transVarLevel >  thresh(1) & transVarLevel <= thresh(2),:) ...
                                 xreg(transVarLevel >  thresh(2),:)}','rls',lower(algo),...
                                {Rest Rest Rest}',{r r r}');
        
        end % switch {checkRest}
        
end % switch {numThresh}

% Residuals
resid = cell(size(res,1),1);
for indRegimes = 1:size(res,1)
	resid{indRegimes} = res(indRegimes).data.ynan - res(indRegimes).data.xnan*res(indRegimes).beta;
end
        
% Compute the log-likelihood
[logL,logCL] = tvarxloglik(resid,commonvar);

% Compute the log of the determinant of the residuals variance-covariance matrix
logDet = tvarxlogdet(resid,commonvar);

% Compute the 
residall = [];
for indRegimes = 1:size(resid,1)
	residall = [residall;resid{indRegimes}]; %#ok<AGROW>
end
ssr = sum(residall(:).^2);
                                            
% Transition variable index
transVarInd = NaN(size(transVarLevel));
switch numThresh
    case 1  % One threshold
        transVarInd(transVarLevel <= thresh) = 1; 
        transVarInd(transVarLevel > thresh)  = 2; 
        
    case 2  % Two thresholds
        transVarInd(transVarLevel <= thresh(1)) = 1; 
        transVarInd(transVarLevel >  thresh(1) & transVarLevel <= thresh(2)) = 2; 
        transVarInd(transVarLevel >  thresh(2)) = 3; 
        
end % switch {numThresh}


% Output into the structure RESULTS
% _________________________________

% Initialize the nested structures
data = ...
    struct('y',     y, ...      % numObs-by-numY, endogenous variables
           'exo',	exo, ...	% numObs-by-numExo, exogenous variables
           'yreg',	[], ...     % (numObs-numLags)-by-numY, regressands (endogenous)
           'xreg',	[], ...     % (numObs-numLags)-by-numX, regressors (both endogenous and exogenous)
           'ynan',	[], ...     % (numObs-numNaN)-by-numY, regressands with NaN removed
           'xnan',	[]);        % (numObs-numNaN)-by-numX, regressors with NaN removed

stats = ...
    struct('r2',	[], ...     % numY-by-1, R-squared
           'rbar2',	[], ...     % numY-by-1, adjusted R-squared
           'vcov',	[], ...     % numX-by-numX, parameters' variance-covariance matrix 
           'sigu',	[], ...     % numY-by-numY, sum of squared residuals    
           'sige',	[]);        % numY-by-numY, residuals variance    

sizes = ...
    struct('numGrid',       numGrid, ...        % scalar, number of points in the grid 
           'numThresh',     numThresh, ...      % scalar, number of threshold
           'numRegimes',    numThresh+1, ...	% scalar, number of regime
           'numLags',       numLags, ...        % scalar, lag length
           'numObs',    	numObs, ...         % scalar, number of observations
           'numObsLag',    	[], ...             % scalar, number of observations adjusted to feed the lags
           'numObsNaNLag',	[], ...             % scalar, number of observations adjusted to feed the lags after NaN are removed
           'numY',          [], ...             % scalar, number of regressands
           'numX',          [], ...             % scalar, total number of regressors
           'numExo',        numExo, ...         % scalar, number of exogenous regressors
           'numNaN',        [], ...             % scalar, number of removed values (NaN)
           'numParams',     [], ...             % scalar, number of parameters (including variance-covariance matrix)
           'numActive',     []);                % scalar, number of active (unrestricted) parameters (including variance-covariance matrix)

optim = ...
    struct('optimcrit',     optimcrit, ...	% string, criterion to optimize
           'gridpt',        gridpt, ...     % numGrid-by-numThresh, matrix of potential threshold
           'logL',          logL, ...       % scalar, log-likelihood evaluated at beta
           'logLloop',      logLloop, ...	% numGrid-by-numThresh, matrix of log-likelihood for each 
           'logDet',        logDet, ...     % scalar, log of the determinant of the residuals variance
           'logDetloop',	logDetloop, ... % numGrid-by-numThresh, matrix of log of the determinant 
           'ssr',           ssr, ...        % scalar, log of the determinant of the residuals variance
           'ssrloop',       ssrloop);   	% numGrid-by-numThresh, matrix of log of the determinant 

threshold = ...
    struct('thresh',        thresh, ...         % numThresh-by-1, thresholds value(s)
           'gamm',          gamm, ...           % vector, prespecified threshold values 
           'threshRatio',   []);                % scalar, proportion of observations in a specific regime

transition = ...
    struct('transVar',      transVar, ...       % logical or double, the input to the function
           'transVarLevel',	transVarLevel, ...	% numObs-by-1, transition variables
           'transVarInd',	transVarInd, ...    % numObs-by-1, vector for the position of regimes 
           'timeDelay',     timeDelay, ...      % scalar, 'time delay' for the threshold variable 
           'trim',          trim);              % scalar, trimming parameter indicating the minimal percentage of observations in each regime 

       
% Fill results structure with field's name and field's value
results(1:numThresh+1,1) = ...
	struct('meth',          'tvarx', ...        % method
           'description',	'Threshold vector autogressive with exogenous variables', ...     % method description
           'fitmeth',       res(1).fitmeth,...	% fitting method
           'algo',          'Grid Search',...	% solution algorithm 
           'optim',         optim, ...          % optimization
           'threshold',     threshold, ...      % threshold
           'transition',    transition, ...     % transition 
           'constant',      cons, ...           % logical, true if constant
           'trend',         trend, ...          % logical, true if linear trend
           'ylog',          ylog, ...           % logical numY-by-1 vector, log the endogenous series before estimating 
           'ydiff',         ydiff, ...          % logical numY-by-1 vector, first difference the endogenous series before estimating 
           'commonvar',     commonvar, ...      % logical, true if all regime   
           'het',           het, ...            % logical, true if residuals are assumed to be hetetoskedastic   
           'corre',         corre, ...          % logical, true if residuals are assumed to be correlated
           'prtr',          prtr, ...           % logical, true if print results
           'graph',         graph, ...          % logical, true if graph fitted values vs observed values
           'vnames',        vnames, ...         % string, vector of variable names (only relevant when prtr = true). If 
           'calstruct',     calstruct, ...      % structure, a calendar structure
           'options',       options, ...        % structure, options provided by the user
           'solmat',        [], ...             % structure, whose fields are matrices of the solution
           'data',          data, ...           % structure, whose fields are matrices of data
           'logCL',         [], ...             % numObs-by-1, vector of log-likelihoods corresponding to beta
           'beta',          [], ...             % structure, estimated coefficients
           'stderr',        [], ...             % structure, standard errors of estimated coefficients
           'tstat',         [], ...             % structure, t-stats of estimated coefficients
           'pval',          [], ...             % structure, p-values of estimated coefficients
           'rest',          [], ...             % structure, restrictions of estimated coefficients
           'yhat',          [], ...             % (numObs-numNaN)-by-numY, fitted values
           'resid',         [], ...             % (numObs-numNaN)-by-numY, residuals
           'stats',         stats, ...          % structure, whose fields are matrices of statistics
           'sizes',         sizes, ...          % structure, whose fields are defined size
           'allNaN',        [], ...             % logical, true if the series is all NaN
           'rankdef',       []);                % logical, true if the series is rank deficient    
           
       
% Fill results structure
% ______________________

% First, fill out results field that are already calculated
for indRegimes = 1:numThresh+1
	
    % Put results inside results structure (help to keep the code readable)
    results(indRegimes).solmat  = res(indRegimes).solmat;
    results(indRegimes).logCL   = logCL{indRegimes};
    results(indRegimes).allNaN	= res(indRegimes).allNaN;
    results(indRegimes).rankdef = res(indRegimes).rankdef;

	results(indRegimes).threshold.threshRatio = res(indRegimes).sizes.numObs/(numObs-numLags);

    results(indRegimes).data.yreg = res(indRegimes).data.y;
    results(indRegimes).data.xreg = res(indRegimes).data.x;
    results(indRegimes).data.ynan = res(indRegimes).data.ynan;
    results(indRegimes).data.xnan = res(indRegimes).data.xnan;

	results(indRegimes).sizes.numObsLag    = res(indRegimes).sizes.numObs;
	results(indRegimes).sizes.numObsNaNLag = res(indRegimes).sizes.numObsNaN;
    results(indRegimes).sizes.numY         = res(indRegimes).sizes.numY;
    results(indRegimes).sizes.numX         = res(indRegimes).sizes.numX;
    results(indRegimes).sizes.numNaN       = res(indRegimes).sizes.numNaN;
    results(indRegimes).sizes.numParams    = res(indRegimes).sizes.numParams + (res(indRegimes).sizes.numY*(res(indRegimes).sizes.numY+1)/2);
    results(indRegimes).sizes.numActive    = res(indRegimes).sizes.numActive + (res(indRegimes).sizes.numY*(res(indRegimes).sizes.numY+1)/2);

end % for {indRegimes}

% Compute useful measures
for indRegimes = 1:numThresh+1
    
    % Compute useful measures 
    results(indRegimes).yhat = results(indRegimes).data.xnan*res(indRegimes).beta;          % fitted values
    results(indRegimes).resid = results(indRegimes).data.ynan - results(indRegimes).yhat;	% residuals

end % for {indRegimes}
    
% Compute overall measures
switch commonvar
    case 0
        % Nothing, everything is regime specific
    case 1
        numObsNaNLag = 0;
        for indRegimes = 1:numThresh+1
            numObsNaNLag = numObsNaNLag + results(indRegimes).sizes.numObsNaNLag;
        end

        yhatall  = NaN(numObsNaNLag,results(1).sizes.numY);
        residall = NaN(numObsNaNLag,results(1).sizes.numY);
        yall     = NaN(numObsNaNLag,results(1).sizes.numY);
        xall     = NaN(numObsNaNLag,results(1).sizes.numX);
        for indRegimes = 1:numThresh+1
            yhatall(transVarInd == indRegimes,:)  = results(indRegimes).yhat;       % fitted values 
            residall(transVarInd == indRegimes,:) = results(indRegimes).resid;      % residuals
            yall(transVarInd == indRegimes,:)     = results(indRegimes).data.ynan;	% regressands
            xall(transVarInd == indRegimes,:)     = results(indRegimes).data.xnan;	% regressors
        end
end

% Loop over all regimes (to get regime specific statistics)
for indRegimes = 1:numThresh+1

    % Compute the variance-covariance matrix
	switch commonvar
        case 0
            % Compute useful measures 
            sigu = results(indRegimes).resid'*results(indRegimes).resid;                                        % residuals sum of square
            sige = results(indRegimes).resid'*results(indRegimes).resid/results(indRegimes).sizes.numObsNaNLag;	% residuals variance
            vcov = varxvcov(results(indRegimes).data.xnan,results(indRegimes).resid,het,corre);
        case 1
            % Compute useful measures 
            sigu = residall'*residall;                  % residuals sum of square
            sige = residall'*residall/numObsNaNLag;     % residuals variance
            vcov = varxvcov(xall,residall,het,corre);
	end
    
    % Compute R-squared
    switch cons
        case 0
            r2 = 1-sum(results(indRegimes).resid.^2)./sum(results(indRegimes).data.ynan.^2);
        case 1    
            ytil = demean(results(indRegimes).data.ynan);
            r2 = 1-sum(results(indRegimes).resid.^2)./sum(ytil.^2);
    end
    r2 = r2';

    % Compute Adjusted R-squared
    switch cons
        case 0
            rbar2 = 1-(sum(results(indRegimes).resid.^2)/(results(indRegimes).sizes.numObs-results(indRegimes).sizes.numX))./...
                (sum(results(indRegimes).data.ynan.^2)/(results(indRegimes).sizes.numObsNaNLag-1));
        case 1    
            ytil  = demean(results(indRegimes).data.ynan);
            rbar2 = 1-(sum(results(indRegimes).resid.^2)/(results(indRegimes).sizes.numObsNaNLag-results(indRegimes).sizes.numX))./...
                (sum(ytil.^2)/(results(indRegimes).sizes.numObsNaNLag-1));
    end
    rbar2 = rbar2';

    % Need to reshape the parameters into a formatted structure. The 
    % parameters are ordered x(-1) y(-1) x(-2) y(-2) ... exo constant trend. 
    tempParam = res(indRegimes).beta';

    % Compute standard errors, t-stats and p-values
    tempStd = reshape(sqrt(diag(vcov)),results(indRegimes).sizes.numX,results(indRegimes).sizes.numY)';
    temptstat = tempParam./tempStd;
    temppval = 2-2*normcdf(abs(temptstat));

    % Preallocate structure for parameters, standard errors, t-stat, p-value and restrictions
    beta = ...
        struct('lags',	{cell(numLags,1)}, ...
               'cons',  [], ...
               'trend', [], ...
               'exo',   []);
    stderr = ...
        struct('lags',	{cell(numLags,1)}, ...
               'cons',  [], ...
               'trend', [], ...
               'exo',   []);
    tstat = ...
        struct('lags',	{cell(numLags,1)}, ...
               'cons',  [], ...
               'trend', [], ...
               'exo',   []);
    pval = ...
        struct('lags',	{cell(numLags,1)}, ...
               'cons',  [], ...
               'trend', [], ...
               'exo',   []);
    rest = ...
        struct('lags',	{cell(numLags,1)}, ...
               'cons',  [], ...
               'trend', [], ...
               'exo',   []);

    % Allocate lags
    for indLags = 1:numLags
        beta.lags{indLags}   = tempParam(:,(indLags-1)*results(indRegimes).sizes.numY+1:indLags*results(indRegimes).sizes.numY);
        stderr.lags{indLags} = tempStd(:,(indLags-1)*results(indRegimes).sizes.numY+1:indLags*results(indRegimes).sizes.numY);
        tstat.lags{indLags}  = temptstat(:,(indLags-1)*results(indRegimes).sizes.numY+1:indLags*results(indRegimes).sizes.numY);
        pval.lags{indLags}   = temppval(:,(indLags-1)*results(indRegimes).sizes.numY+1:indLags*results(indRegimes).sizes.numY);
        rest.lags{indLags}   = ARsolve(:,:,indLags);
    end

    % Allocate constant
    switch cons
        case 0
            % No constant
        case 1
            beta.cons   = tempParam(:,end-double(trend));
            stderr.cons = tempStd(:,end-double(trend));
            tstat.cons  = temptstat(:,end-double(trend));
            pval.cons   = temppval(:,end-double(trend));
            rest.cons   = Csolve;
    end

    % Allocate linear trend
    switch trend
        case 0
            % No linear trend
        case 1
            beta.trend   = tempParam(:,end);
            stderr.trend = tempStd(:,end);
            tstat.trend  = temptstat(:,end);
            pval.trend   = temppval(:,end);         
            rest.trend   = Tsolve;         
    end

    % Allocate exogenous regressors
    checkExo = (size(tempParam,2) ~= double(cons) + double(trend) + numLags*results(indRegimes).sizes.numY);
    switch checkExo
        case 0
            % No exogenous regressors
        case 1
            beta.exo   = tempParam(:,numLags*results(indRegimes).sizes.numY+1:end-double(cons)-double(trend));
            stderr.exo = tempStd(:,numLags*results(indRegimes).sizes.numY+1:end-double(cons)-double(trend));
            tstat.exo  = temptstat(:,numLags*results(indRegimes).sizes.numY+1:end-double(cons)-double(trend));
            pval.exo   = temppval(:,numLags*results(indRegimes).sizes.numY+1:end-double(cons)-double(trend));
            rest.exo   = EXOsolve;
    end

    % Put results inside results structure (help to keep the code readable)
    results(indRegimes).beta   = beta;
    results(indRegimes).stderr = stderr;	
    results(indRegimes).tstat  = tstat;
    results(indRegimes).pval   = pval;
    results(indRegimes).rest   = rest;

	results(indRegimes).stats.sigu	= sigu;
    results(indRegimes).stats.sige	= sige;
    results(indRegimes).stats.r2    = r2;
    results(indRegimes).stats.rbar2	= rbar2;
    results(indRegimes).stats.vcov	= vcov;

end % for {indRegimes}


% Print results      
% _____________

switch prtr
    case 0
        % No print
    case 1
        switch ~isempty(vnames)
        	case 0
                tvarxprt(results);
            case 1
                tvarxprt(results,vnames)
        end % switch {~isempty(vnames)}
end % switch {prtr}
    

% Graph fitted vs observed values      
% _______________________________

switch graph
    case 0
        % No graph
    case 1
        switch ~isempty(vnames)
            case 0
                tvarxplt(results);
            case 1
                tvarxplt(results,vnames)
        end % switch {~isempty(vnames)}        
end % switch {graph}


end % function {tvarx}


% ---------------------------------
function [y,numLags,numThresh,transVar,timeDelay] = checkInput(y,numLags,numThresh,transVar,timeDelay)
% checkInput: Local function to check the validity of required inputs

% y
validateattributes(y,{'numeric'},{'2d','real'},'tvarx','y',1);

% numLags
validateattributes(numLags,{'numeric'},{'scalar','integer','positive','finite'},'tvarx','numLags',2);

% numThresh
validateattributes(numThresh,{'numeric'},{'scalar','integer','positive','>=',1,'<=',2},'tvarx','numThresh',3);   

% transVar
validateattributes(transVar,{'logical','numeric'},{},'tvarx','transVar',4);
if islogical(transVar)                                                                                              
    validateattributes(transVar,{'logical'},{'numel',size(y,2)},'tvarx','transVar',4);
    if sum(double(transVar)) > 1
        error('tvarx:InputError','tvarx: currently, for transVar, only a single value is implemented')
    elseif sum(double(transVar)) == 0
        error('tvarx:InputError','tvarx: at least one transition variable must be declared')
    end
elseif isnumeric(transVar)
    validateattributes(transVar,{'numeric'},{'real','nrows',size(y,1),'ncols',1},'tvarx','transVar',4);
end

% timeDelay
validateattributes(timeDelay,{'numeric'},{'real','scalar','integer','nonnegative'},'tvarx','timeDelay',5);   
if islogical(transVar)                                                                                              
    validateattributes(timeDelay,{'numeric'},{'>=',1,'<=',numLags},'tvarx','timeDelay',5);
elseif isnumeric(transVar)
    validateattributes(timeDelay,{'numeric'},{'>=',0,'<=',numLags},'tvarx','timeDelay',5);
end


end % subfunction {checkInput}


% ---------------------------------
function [cons,trend,exo,ylog,ydiff,Rest,algo,commonvar,het,corre,optimcrit,numGrid,gamm,trim,around,prtr,graph,vnames,calstruct] = checkOptions(options,y,numThresh) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Check if exo is a user-provided options
switch isfield(options,'exo')
    case 1
        numExo = size(options.exo,2);
    case 0
        numExo = 0;
end

% Get default or user-provided options
getOptions(options, ...
	'cons',         false, ...                          % constant
    'trend',        false, ...                          % linear trend
    'exo',          double.empty(size(y,1),0), ...      % exogenous regressors
    'ylog',         false(1,size(y,2)), ...             % log the endogenous variables
    'ydiff',        false(1,size(y,2)), ...             % first difference the endogenous variables
    'Rest',         struct.empty, ...                   % Restrictions structure
    'algo',         'qr', ...                           % way to compute the MLE estimator of beta
    'commonvar',	false, ...                          % share the same variance-covariance matrix
    'het',          false, ...                          % heteroskedasticity
    'corre',        false, ...                          % correlated residuals
    'optimcrit',    'loglik', ...                       % optimization criterion
    'numGrid',      50, ...                             % number of equally spaced elements of the grid 
    'gamm',         double.empty, ...                   % prespecified threshold values
    'trim',         15, ...                             % trimming parameter indicating the minimal percentage of observations in each regime
    'around',       false, ...                          % perform a second optimization round 
    'prtr',         false, ...                          % print results
	'graph',        false, ...                          % graph fitted vs observed values
    'vnames',       char.empty(size(y,2)+numExo,0), ...	% vector of variable names (only relevant when print = true)
    'calstruct',	struct.empty);                      % a calendar structure

% cons
validateattributes(cons,{'logical'},{'numel',1},'tvarx','options.cons',6);	

% trend
validateattributes(trend,{'logical'},{'numel',1},'tvarx','options.trend',6);

% ylog
validateattributes(ylog,{'logical'},{'numel',size(y,2)},'tvarx','options.ylog',6);	%#ok<NODEF>
switch any(ylog)
    case 1
        warning('tvarx:InputProblem','tvarx: Do not use "ylog" for the moment. The other functions are not yet updated to reflect these new options. "ylog" will be set to false(1,size(y,2)).')
        ylog = false(1,size(y,2));
    case 0
        % Nothing
end

% ydiff
validateattributes(ydiff,{'logical'},{'numel',size(y,2)},'tvarx','options.ydiff',6); %#ok<NODEF>
switch any(ydiff)
    case 1
        warning('tvarx:InputProblem','tvarx: Do not use "ydiff" for the moment. The other functions are not yet updated to reflect these new options. "ydiff" will be set to false(1,size(y,2))')
        ydiff = false(1,size(y,2));
    case 0
        % Nothing
end


% exo
validateattributes(exo,{'numeric'},{'2d','real','nrows',size(y,1)},'tvarx','options.exo',6);  

% Rest
validateattributes(Rest,{'struct'},{},'tvarx','options.Rest',6);

% algo
validateattributes(algo,{'char'},{},'varx','options.algo',3);
validatestring(lower(algo),{'explicit' 'qr' 'svd'},'tvarx','options.algo',6);

% commonvar
validateattributes(commonvar,{'logical'},{'numel',1},'tvarx','options.commonvar',6);

% het
validateattributes(het,{'logical'},{'numel',1},'tvarx','options.het',6);

% corre
validateattributes(corre,{'logical'},{'numel',1},'tvarx','options.corre',6);

% optimcrit
validateattributes(optimcrit,{'char'},{},'tvarx','options.optimcrit',6);
validatestring(lower(optimcrit),{'loglik' 'logdet' 'ssr'},'tvarx','options.optimcrit',6);

% numGrid
validateattributes(numGrid,{'numeric'},{'real','scalar','integer','positive'},'tvarx','options.numGrid',6);

% gamm
if ~isempty(gamm)  
    validateattributes(gamm,{'numeric'},{'real','numel',numThresh},'tvarx','options.gamm',6);
end

% trim
validateattributes(trim,{'numeric'},{'real','positive','scalar','>=',10,'<=',50},'tvarx','options.trim',6); 

% around
validateattributes(around,{'logical'},{'numel',1},'tvarx','options.around',6);	%#ok<NODEF>
switch around
    case 1
        warning('tvarx:InputProblem','tvarx: The options "around" is not working very well. "around" will be set to false. Increase the number of grid points instead.')
        around = false;
    case 0
        % Nothing
end

% prtr
validateattributes(prtr,{'logical'},{'numel',1},'tvarx','options.prtr',6);

% graph
validateattributes(graph,{'logical'},{'numel',1},'tvarx','options.graph',6);

% vnames
validateattributes(vnames,{'char'},{'nrows',size(y,2)+size(exo,2)},'tvarx','options.vnames',6);

% calstruct
validateattributes(calstruct,{'struct'},{},'tvarx','options.calstruct',6);


end % subfunction {checkOptions}

