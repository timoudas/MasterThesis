function results = varx(y,numLags,options)
% VARX: Estimate a unrestricted VARX model by ordinary least squares (equivalent to MLE) or a 
%		restricted VARX model by generalized least squares.
%
% Syntax:
%
%       results = varx(y,numLags)
%       results = varx(y,numLags,options)
%       results = varx(y,numLags,struct('cons',true',trend',true))
%
% Description:
%
%       Estimate a vector autogressive with (potentially) exogenous variables model by least squares of the following form:
%
%       y(t) = A1*y(t-1) + ... + Ap*y(t-p) + B*x(t) + C*D(t) + e(t)
%
%		The unrestricted VARX model is estimated by ordinary least squares (NB: the ordinary least squares
%       is equivalent to maximum likelihood when one iteration is computed.) and the  
%		restricted VARX model is estimated by generalized least squares.
%
%       VARX treats NaNs in 'y' and 'exo' as missing values, and removes them.
%
%       VARX can automatically include a constant and a linear trend. Other deterministic regressors must be included into 
%       the matrix of exogenous regressors 'exo'.
%
% Input Arguments:
%
%       y               -	[double]	numObs-by-numY, dependent variable matrix (or response observations, or 
%                                       regressands)
%       numLags         -   [integer]   scalar, lag length
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       Model Specification options:
%       	cons        -   [logical]   a constant
%                                           true	include a constant (default)
%       	                                false	does not include a constant
%       	trend       -   [logical]   a linear trend
%       	                                true	include a linear trend 
%       	                                false	does not include a linear trend (default)
%       	exo         -	[double]	numObs-by-numVars, independent variables matrix (regressors). Its a design 
%       	                            matrix, with rows corresponding to observations and columns to predictor 
%       	                            variables graph 
%       	                            	Default: empty
%       	ylog        -   [logical]   a logical numY-by-1 vector, log the endogenous series before estimating 
%                                       the VARX: true means apply the log, false means no log is taken 
%       	                            	Default: No log (logical(zeros(numY,1)))	
%       	ydiff       -   [logical]   a logical numY-by-1 vector, first difference the endogenous series before estimating 
%                                       the VARX: true means apply the difference, false means no difference is taken 
%       	                            	Default: No difference (logical(zeros(numY,1)))	
%			Rest		-	[struct]	structure of numerical restrictions. All restrictions are defined this way:
%											1 - NaN means no restrictions applied on this parameters
%											2 - a real numeric value to applied this numerical restriction on the specific parameter
%                                       The structure has these specifics fields:
%       	                            	ARsolve     a numY-by-numY-by-numLags matrix (default: NaN(numY,numY,numLags)),
%                                                       if more than numLags layer are provided, only the first numLags will be kept
%											Csolve      a numY-by-1 matrix (default: if cons = true, a NaN(numY,1))
%											Tsolve      a numY-by-1 matrix (default: if trend = true, a NaN(numY,1))
%											EXOsolve	a numY-by-numX matrix (default: if exo not empty, a NaN(numY,numX))
%           algo        -   [char]      string, way to compute the MLE estimator of beta:
%                                           'explicit'  explicit formula for global minimum
%                                           'qr'        QR decomposition (default)
%                                           'svd'       SVD decomposition
%       Variance-Covariance options:
%       	het         -   [logical]   type of covariance estimator
%       	                                true    heteroskedastic
%       	                                false	homoskedastic (default)
%       	corre       -   [logical]   the assumed structure of the error covariance matrix
%       	                                true    correlated residuals
%       	                                false	uncorrelated residuals (default)
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
%           calstruct	-   [struct]    a calendar struct, as returned by cal.m
%
% Output Arguments:
%
%       results, a structure with the following fields:         
%           meth        -   [char]      method
%           description	-   [char]      method description
%           fitmeth     -   [char]      fitting method
%           algo        -   [char]      string, solution algorithm
%           solmat      -   [struct]    structure, whose fields are matrices of the solution
%                                       If algo = 'explicit':
%                                       	nothing
%                                       If algo = 'qr':
%                                           Q           - 	[double]    (numObs-numNaN)-by-numX, Q from the QR decomposition of the design matrix
%                                           R           - 	[double]    numX-by-numX, R from the QR Decomposition of the design matrix    
%                                       	perm        -   [double]    1-by-numX, permutation information
%                                       If algo = 'svd':
%                                           U           - 	[double]    (numObs-numNaN)-by-(numObs-numNaN), U from the SVD decomposition
%                                                                       of the design matrix    
%                                           V           - 	[double]    numX-by-numX, V from the SVD decomposition of the design matrix
%                                           S           -   [double]    (numObs-numNaN)-by-numX, rectangular diagonal matrix with 
%                                                                       nonnegative real numbers on the diagonal
%           logL        -   [double]    scalar, log-likelihood evaluated at beta
%           logCL       -   [double]    numObs-by-1, vector of log-likelihoods corresponding to beta
%           data        -	[struct]	structure, whose fields are matrices are the data
%                                       	y           -   [double]    numObs-by-numY, endogenous variables
%                                           exo         -   [double]    numObs-by-numExo, exogenous variables
%                                           yreg        -   [double]    (numObs-numLags)-by-numY, regressands (endogenous)
%                                           xreg        -   [double]    (numObs-numLags)-by-numX, regressors (both endogenous and exogenous)
%                                           ynan        -   [double]    (numObs-numNaN)-by-numY, regressands with NaN removed
%                                           xnan        -   [double]    (numObs-numNaN)-by-numX, regressors with NaN removed
%           beta        -   [struct]    structure, estimated coefficients 
%           stderr      -   [struct]    structure, standard errors of estimated coefficients
%           tstat       -   [struct]    structure, t-stats of estimated coefficients
%           pval        -   [struct]    structure, p-values of estimated coefficients
%           rest        -   [struct]    structure, restrictions on estimated coefficients
%           yhat        -   [double]    (numObs-numNaN)-by-numY, fitted values
%           resid       -   [double]    (numObs-numNaN)-by-numY, residuals
%           stats       -   [struct]    structure, whose fields are matrices are the statistics
%                                       	r2          -   [double]    numY-by-1, R-squared
%                                           rbar2       -   [double]    numY-by-1, adjusted R-squared
%                                           vcov        -   [double]    numX-by-numX, parameters' variance-covariance matrix 
%                                           sigu        -   [double]    numY-by-numY, sum of squared residuals    
%                                           sige        -   [double]    numY-by-numY, residuals variance    
%           sizes       -   [struct]    structure, whose fields are defined size
%                                       	numObs      -   [double]    scalar, number of observations
%                                           numLags     -   [double]    scalar, lag length
%                                           numObsLag   -   [double]    scalar, number of observations adjusted to feed the lags
%                                           numObsNaNLag-   [double]    scalar, number of observations adjusted to feed the lags after NaN are removed
%                                           numY        -   [double]    scalar, number of regressands
%                                           numX        -   [double]    scalar, total number of regressors
%                                           numExo      -   [double]    scalar, number of exogenous regressors
%                                           numNaN      -   [double]    scalar, number of removed values (NaN)
%                                           numParams   -   [double]    scalar, number of parameter
%                                           numActive   -   [double]    scalar, number of unrestricted parameter
%           constant	-   [logical]   true if constant
%           trend       -   [logical]   true if linear trend
%           ylog        -   [logical]	logical numY-by-1 vector, log the endogenous series before estimating 
%           ydiff       -   [logical]	logical numY-by-1 vector, first difference the endogenous series before estimating 
%           het         -   [logical]   true is residuals are assumed to be hetetoskedastic   
%           corre       -   [logical]   true is residuals are assumed to be correlated
%           prtr        -   [logical]   true if print results
%           graph       -   [logical]   true if graph fitted values vs observed values
%           vnames      -   [char]      string, vector of variable names (only relevant when prtr = true). 
%           calstruct	-   [struct]    structure, a calendar structure
%           allNaN      -   [logical]   true if the series is all NaN
%           rankdef     -   [logical]   true if the series is rank deficient 
%           options     -   [struct]    structure, user-provided options
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, lagmat, mvlinreg, demean, varxvcov, varxloglik, varxprt, varxplt
%       
%
% References:
%
%       none
%
% Notes:
%
%       A difference in the estimation of the parameters with other VARX package when lag is > 2 is possible. 
%       This is because I initialize the trend with the number of lags (i.e., when numLags = 2, the trend is [2 3 ...numObs]), 
%       while other package may initialize the trend always with 1. 
%
%       Do not use ylog and ydiff for the moment. The other functions are not yet updated to reflect these new options.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014
%           based on vare.m from James P. LeSage's Econometrics Toolbox:
%               James P. LeSage, Dept of Economics
%               University of Toledo
%               2801 W. Bancroft St,
%               Toledo, OH 43606
%               jlesage@spatial-econometrics.com
%           and vectorar.m from Kevin Sheppard MFE Toolbox:
%               kevin.sheppard@economics.ox.ac.uk


% Input and Output arguments checking
% ___________________________________

narginchk(2,3); nargoutchk(1,1);

[y,numLags] = checkInput(y,numLags);


% Options
% _______

if nargin < 3
    options = struct.empty;
elseif nargin == 3
	validateattributes(options,{'struct'},{},'varx','options',3);
end

[cons,trend,exo,ylog,ydiff,Rest,algo,het,corre,prtr,graph,vnames,calstruct] = checkOptions(options,y);


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
    validateattributes(ARsolve,{'numeric'},{'size',[size(y,2) size(y,2) numLags]},'varx','options.Rest.ARsolve',3);
elseif size(ARsolve,3) > numLags
    ARsolve = ARsolve(:,:,1:numLags);
    validateattributes(ARsolve,{'numeric'},{'size',[size(y,2) size(y,2) numLags]},'varx','options.Rest.ARsolve',3);
end

% Csolve
validateattributes(Csolve,{'numeric'},{'size',[size(y,2) double(cons)]},'varx','options.Rest.Csolve',3);

% Tsolve
validateattributes(Tsolve,{'numeric'},{'size',[size(y,2) double(trend)]},'varx','options.Rest.Tsolve',3);

% EXOsolve
validateattributes(EXOsolve,{'numeric'},{'size',[size(y,2) size(exo,2)]},'varx','options.Rest.EXOsolve',3);


% Check inputs, form inputs and restrictions matrices 
% ___________________________________________________

% Check endogenous inputs
numObs = size(y,1); % Define numObs, the number of observations
Restrictions = [Restrictions reshape(ARsolve,[size(y,2) size(y,2)*numLags 1])];

% Check exogenous inputs
switch ~isempty(exo)
    case 0
        % No exogenous regressors, than no EXO restrictions
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
switch all(all(isnan(Restrictions)))
	case 1	% No restrictions
        % Nothing
    case 0  % Restrictions
		Rest = eye(size(yreg,2)*size(xreg,2));
		r = Restrictions(:);

		Rest = Rest(any(isnan(r),2),:).';
		r(isnan(r)) = 0;
end


% Estimation
% __________

% Convert yreg and xreg to cell array, because mvlinreg only accept cells
switch all(all(isnan(Restrictions)))
	case 1	% No restrictions, use OLS
		res = mvlinreg({yreg},{xreg},'ols',algo);
    case 0	% Restrictions, use RLS
		res = mvlinreg({yreg},{xreg},'rls',algo,{Rest},{r});
end


% Compute measures and statistics 
% _______________________________

yhat  = res.data.xnan*res.beta;                 % fitted values 
resid = res.data.ynan - yhat;                   % residuals
sigu  = resid'*resid;                           % residuals sum of square
sige  = resid'*resid/res.sizes.numObsNaN;       % residuals variance

% Compute the log-likelihood and individual log-likelihoods
[logL,logCL] = varxloglik(resid);

% Compute the variance-covariance matrix of estimated coefficients
vcov = varxvcov(res.data.xnan,resid,het,corre);

% Compute R-squared
switch cons
    case 0
        r2 = 1-sum(resid.^2)./sum(res.data.ynan.^2);
    case 1    
        ytil = demean(res.data.ynan);
        r2 = 1-sum(resid.^2)./sum(ytil.^2);
end
r2 = r2';

% Compute Adjusted R-squared
switch cons
    case 0
        rbar2 = 1-(sum(resid.^2)/(res.sizes.numObs-res.sizes.numX))./(sum(res.data.ynan.^2)/(res.sizes.numObs-1));
    case 1    
        ytil = demean(res.data.ynan);
        rbar2 = 1-(sum(resid.^2)/(res.sizes.numObs-res.sizes.numX))./(sum(ytil.^2)/(res.sizes.numObs-1));
end
rbar2 = rbar2';


% Compute standard errors, t-stats and p-values 
% _____________________________________________

% Need to reshape the parameters into a formatted structure. The 
% parameters are ordered x(-1) y(-1) x(-2) y(-2) ... exo constant trend. 
tempParam = res.beta.';

% Compute standard errors, t-stats and p-values
tempStd   = reshape(sqrt(diag(vcov)),res.sizes.numX,res.sizes.numY)';
temptstat = tempParam./tempStd;
temppval  = 2-2*normcdf(abs(temptstat));


% Organize parameters, standard errors, t-stat, p-value and restrictions into matrices and structure
% __________________________________________________________________________________________________

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
	beta.lags{indLags}   = tempParam(:,(indLags-1)*res.sizes.numY+1:indLags*res.sizes.numY);
    stderr.lags{indLags} = tempStd(:,(indLags-1)*res.sizes.numY+1:indLags*res.sizes.numY);
    tstat.lags{indLags}  = temptstat(:,(indLags-1)*res.sizes.numY+1:indLags*res.sizes.numY);
    pval.lags{indLags}   = temppval(:,(indLags-1)*res.sizes.numY+1:indLags*res.sizes.numY);
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
checkExo = (size(tempParam,2) ~= double(cons) + double(trend) + numLags*res.sizes.numY);
switch checkExo
	case 0
        % No exogenous regressors
	case 1
        beta.exo   = tempParam(:,numLags*res.sizes.numY+1:end-double(cons)-double(trend));
        stderr.exo = tempStd(:,numLags*res.sizes.numY+1:end-double(cons)-double(trend));
        tstat.exo  = temptstat(:,numLags*res.sizes.numY+1:end-double(cons)-double(trend));
        pval.exo   = temppval(:,numLags*res.sizes.numY+1:end-double(cons)-double(trend));
        rest.exo   = EXOsolve;
end


% Output into the structure RESULTS
% _________________________________

% Initialize the nested structures
data = ...
    struct('y',     y, ...              % numObs-by-numY, endogenous variables
           'exo',	exo, ...            % numObs-by-numExo, exogenous variables
           'yreg',	res.data.y, ...     % (numObs-numLags)-by-numY, regressands (endogenous)
           'xreg',	res.data.x, ...     % (numObs-numLags)-by-numX, regressors (both endogenous and exogenous)
           'ynan',	res.data.ynan, ...	% (numObs-numNaN)-by-numY, regressands with NaN removed
           'xnan',	res.data.xnan);     % (numObs-numNaN)-by-numX, regressors with NaN removed

stats = ...
    struct('r2',	r2, ...     % numY-by-1, R-squared
           'rbar2',	rbar2, ...	% numY-by-1, adjusted R-squared
           'vcov',	vcov, ...	% numX-by-numX, parameters' variance-covariance matrix 
           'sigu',	sigu, ...	% numY-by-numY, sum of squared residuals    
           'sige',	sige);      % numY-by-numY, residuals variance    

sizes = ...
    struct('numObs',    	numObs, ...                 % scalar, number of observations
           'numLags',       numLags, ...                % scalar, lag length
           'numObsLag',    	res.sizes.numObs, ...       % scalar, number of observations adjusted to feed the lags
           'numObsNaNLag',	res.sizes.numObsNaN, ...	% scalar, number of observations adjusted to feed the lags after NaN are removed
           'numY',          res.sizes.numY, ...         % scalar, number of regressands
           'numX',          res.sizes.numX, ...         % scalar, total number of regressors
           'numExo',        numExo, ...                 % scalar, number of exogenous regressors
           'numNaN',        res.sizes.numNaN, ...       % scalar, number of removed values (NaN)
           'numParams',     res.sizes.numParams+(res.sizes.numY*(res.sizes.numY+1)/2), ...	% scalar, number of parameter (including variance-covariance matrix)            %
           'numActive',     res.sizes.numActive+(res.sizes.numY*(res.sizes.numY+1)/2));     % scalar, number of unrestricted parameter (including variance-covariance matrix)	%

% Fill results structure with field's name and field's value
results = ...
	struct('meth',          'varx', ...         % method
           'description',	'Vector autogressive with exogenous variables', ...     % method description
           'fitmeth',       res.fitmeth,...     % fitting method
           'algo',          res.algo,...        % solution algorithm          
           'solmat',        res.solmat, ...     % structure, whose fields are matrices of the solution
           'logL',          logL, ...           % scalar, log-likelihood evaluated at beta
           'logCL',         logCL, ...          % numObs-by-1, vector of log-likelihoods corresponding to beta
           'data',          data, ...           % structure, whose fields are matrices of data
           'beta',          beta, ...           % structure, estimated coefficients
           'stderr',        stderr, ...         % structure, standard errors of estimated coefficients
           'tstat',         tstat, ...          % structure, t-stats of estimated coefficients
           'pval',          pval, ...           % structure, p-values of estimated coefficients
           'rest',          rest, ...           % structure, restrictions of estimated coefficients
           'yhat',          yhat, ...           % (numObs-numNaN)-by-numY, fitted values
           'resid',         resid, ...          % (numObs-numNaN)-by-numY, residuals
           'stats',         stats, ...          % structure, whose fields are matrices of statistics
           'sizes',         sizes, ...          % structure, whose fields are defined size
           'constant',      cons, ...           % logical, true if constant
           'trend',         trend, ...          % logical, true if linear trend
           'ylog',          ylog, ...           % logical numY-by-1 vector, log the endogenous series before estimating 
           'ydiff',         ydiff, ...          % logical numY-by-1 vector, first difference the endogenous series before estimating 
           'het',           het, ...            % logical, true if residuals are assumed to be hetetoskedastic   
           'corre',         corre, ...          % logical, true if residuals are assumed to be correlated
           'prtr',          prtr, ...           % logical, true if print results
           'graph',         graph, ...          % logical, true if graph fitted values vs observed values
           'vnames',        vnames, ...         % string, vector of variable names (only relevant when prtr = true). If 
           'calstruct',     calstruct, ...      % structure, a calendar structure
           'allNaN',        res.allNaN, ...     % logical, true if the series is all NaN
           'rankdef',       res.rankdef, ...    % logical, true if the series is rank deficient    
           'options',       options);           % structure, options provided by the user

       
% Print results      
% _____________

switch prtr
    case 0
        % No print
    case 1
        %checkVnames = ~isempty(vnames);
            switch ~isempty(vnames)
                case 0
                    varxprt(results);
                case 1
                    varxprt(results,vnames)
            end % switch {checkVnames}
end % switch {prtr}
    

% Graph fitted vs observed values      
% _______________________________

switch graph
    case 0
        % No graph
    case 1
        %checkVnames = ~isempty(vnames);
            switch ~isempty(vnames)
                case 0
                    varxplt(results);
                case 1
                    varxplt(results,vnames)
            end % switch {checkVnames}        
end % switch {graph}


end % function {varx}


% ---------------------------------
function [y,numLags] = checkInput(y,numLags)
% checkInput: Local function to check the validity of required inputs

% y
validateattributes(y,{'numeric'},{'2d','real'},'varx','y',1);

% numLags
validateattributes(numLags,{'numeric'},{'scalar','integer','positive','finite'},'varx','numLags',2);


end % subfunction {checkInput}


% ---------------------------------
function [cons,trend,exo,ylog,ydiff,Rest,algo,het,corre,prtr,graph,vnames,calstruct] = checkOptions(options,y) %#ok<STOUT>
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
	'cons',         true, ...                           % constant
    'trend',        false, ...                          % linear trend
    'exo',          double.empty(size(y,1),0), ...      % exogenous regressors
    'ylog',         false(1,size(y,2)), ...             % log the endogenous variables
    'ydiff',        false(1,size(y,2)), ...             % first difference the endogenous variables
    'Rest',         struct.empty, ...                   % restrictions structure
    'algo',         'qr', ...                           % way to compute the MLE estimator of beta
    'het',          false, ...                          % heteroskedasticity
    'corre',        false, ...                          % correlated residuals
    'prtr',         false, ...                          % print results
	'graph',        false, ...                          % graph fitted vs observed values
    'vnames',       char.empty(size(y,2)+numExo,0), ...	% vector of variable names (only relevant when print = true)
    'calstruct',	struct.empty);                      % a calendar structure

% cons
validateattributes(cons,{'logical'},{'numel',1},'varx','options.cons',3);	

% trend
validateattributes(trend,{'logical'},{'numel',1},'varx','options.trend',3);

% ylog
validateattributes(ylog,{'logical'},{'size',[1 size(y,2)]},'varx','options.ylog',3);

% ydiff
validateattributes(ydiff,{'logical'},{'size',[1 size(y,2)]},'varx','options.ydiff',3);

% exo
validateattributes(exo,{'numeric'},{'2d','real','nrows',size(y,1)},'varx','options.exo',3);  

% Rest
validateattributes(Rest,{'struct'},{},'varx','options.Rest',3);

% algo
validateattributes(algo,{'char'},{},'varx','options.algo',3);
validatestring(algo,{'explicit' 'qr' 'svd'},'varx','options.algo',3);

% het
validateattributes(het,{'logical'},{'numel',1},'varx','options.het',3);

% corre
validateattributes(corre,{'logical'},{'numel',1},'varx','options.corre',3);

% prtr
validateattributes(prtr,{'logical'},{'numel',1},'varx','options.prtr',3);

% graph
validateattributes(graph,{'logical'},{'numel',1},'varx','options.graph',3);

% vnames
validateattributes(vnames,{'char'},{'nrows',size(y,2)+size(exo,2)},'varx','options.vnames',3);

% calstruct
validateattributes(calstruct,{'struct'},{},'varx','options.calstruct',3);


end % subfunction {checkOptions}

