function [icVar,infocrit,tresh] = tvarxic(y,nARmax,nThreshmax,transVar,options_ic,options_tvarx)
% TVARXIC: Evaluates criteria for selecting the order, time delay and threshold number of a TVARX model.
%
% Syntax:
%
%       ic = tvarxic(y,nARmax,nThreshMax,transVar,options_ic,options_tvarx)
%
% Description:
%
%       Evaluates criteria for selecting the order of a TVARX model. The information criteria will determine the optimal lag length, the
%       optimal number of threshold and the optimal time delay.
%
%       TVARXIC treats NaNs in 'y' and 'exo' as missing values, and removes them.
%
%       TVARXIC can automatically include a constant and a linear trend. Other deterministic regressors must be included into 
%       the matrix of exogenous regressors 'exo'.
%
% Input Arguments:
%
%       y           -	[double]            numObs-by-numY, dependent variable matrix (or response observations, or regressands)
%       nARmax      -   [integer]           scalar, maximum lag length
%       nThreshMax	-   [integer]           scalar, maximum number of thresholds. Currently, the number of thresholds could be 1 (= 2 regimes)
%                                           or 2 (= 3 regimes).
%       transVar	-   [logical,double]    transVar is either:
%                                               1 - a logical numY-by-1 vector indicating which variables (a single 
%                                               	value or a combination of endogenous and exogenous variables with same lag 
%                                                   order) to take for the transition variables.
%                                                   ex: if y is a 500-by-4 , if transVar is logical([1 0 0 0]), 
%                                                   	then the transition variable is the first one in y. 
%                                               2 - a numObs-by-1 vector of exogenous (included in exo) or external (not included in y or exo) 
%                                               	transition variable.
%
% Optional Input Arguments:
%
%       options_ic, a structure for the information criterion computation with the following possible fields and possible values:
%       infoc       -   [char]          a string
%                                           AIC     for Akaike (default)
%                                           AICc	for Corrected Akaike
%                                           SBIC	for Schwarz-Bayesian
%                                           HQIC	for Hannan-Quinn
%       graph       -	[logical]       print information
%                                           true	graph and print information (default)
%                                           false	do not graph and print information
%       options_tvarx, a structure for the estimation of the TVARX process with the following possible fields and possible values:
%       All the options fields and values in varx.m (see tvarx.m for details)
%
% Output Arguments:
%
%       icVar       -   [scalar]    scalar, optimal lag length according to the chosen criteria
%       infocrit	-   [struct]    structure, all information criteria computed for all lag length, all time delay and all threshold number
%       tresh       -   [cell]      cell array, all threshold computed for all lag length, all time delay and all threshold number
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, tvarx, infocriterion
%
% References:
%
%       none
%
% Notes:
%
%       The information criteria computed with this program (AIC, AICc, SBIC and HQIC) are not convergent for a TVARX process, since they 
%       all need a continuous likelihood function, and the likelihood function of a TVARX is discontinuous. However, some Monte Carlo 
%       studies show that those information criterai still give good results, even in small sample. Therefore, they should be interpreted
%       as a "good" indication, but never as a estimated value that converged to its true value.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(4,6); nargoutchk(1,3);

[y,nARmax,nThreshmax,transVar] = checkInput(y,nARmax,nThreshmax,transVar);


% Options for information criterion
% _________________________________

if nargin < 5
    options_ic = struct.empty;
elseif nargin >= 5
	validateattributes(options_ic,{'struct'},{},'tvarxic','options_ic',5);
end

[infoc,graph] = checkOptions(options_ic);


% Options for TVARX estimation
% ___________________________

if nargin < 6
    options_tvarx = struct.empty;
elseif nargin == 6
	validateattributes(options_tvarx,{'struct'},{},'varxic','options_tvarx',6);
end

% Set options by default
checkOptim = isfield(options_tvarx,'optimcrit');
switch checkOptim
	case 0
        options_tvarx.optimcrit = 'loglik';
	case 1
        % Nothing
end % switch {checkExo}


% Check provided exogenous inputs, if any
% _______________________________________

isExo = isfield(options_tvarx,'exo');
switch isExo
	case 0
    	% No exogenous variables
	case 1
        validateattributes(options_tvarx.exo,{'numeric'},{'2d','real','nrows',size(y,1)},'tvarxic','options_varx.exo',6);  
        exo = options_tvarx.exo;
end % switch {checkExo}


% Preallocate arrays
% __________________

isNum = double(isnumeric(transVar));

data      = cell(nARmax,nARmax+isNum,nThreshmax);
numParams = zeros(nARmax,nARmax+isNum,nThreshmax);
numObs    = zeros(nARmax,nARmax+isNum,nThreshmax);
tresh     = cell(nARmax,nARmax+isNum,nThreshmax);


% Loop over models
% ________________

for indAR = 1:nARmax
    for indTime = 1-isNum:indAR
        for indTresh = 1:nThreshmax
            switch isExo
                case 0  % No exogenous variables
                    % Nothing
                case 1  % Exogenous variables
                    options_tvarx.exo = exo(nARmax+1-indAR:end,:);
            end % switch {checkExo}
            
            % Estimate the model
            if islogical(transVar)
                results = tvarx(y(nARmax+1-indAR:end,:),indAR,indTresh,transVar,indTime,options_tvarx);
            elseif isnumeric(transVar)
                results = tvarx(y(nARmax+1-indAR:end,:),indAR,indTresh,transVar(nARmax+1-indAR:end,:),indTime,options_tvarx);
            end
            
            % Extract results
            if strcmp(options_tvarx.optimcrit,'loglik')
                data{indAR,indTime+isNum,indTresh} = results(1).optim.logL;
            elseif strcmp(options_tvarx.optimcrit,'logdet')
                tempresid = [];
                for indRegime = 1:size(results,1)
                    tempresid = [tempresid;results(indRegime).resid]; %#ok<AGROW>
                end
                data{indAR,indTime+isNum,indTresh} = permute(tempresid,[3 1 2]);
            end
            
            for indRegime = 1:size(results,1)
            	numParams(indAR,indTime+isNum,indTresh) = numParams(indAR,indTime+isNum,indTresh) + results(indRegime).sizes.numActive;
            end
                
            numObs(indAR,indTime+isNum,indTresh) = results(1).sizes.numObsNaNLag;
            
            % Extract threshold
            tresh{indAR,indTime+isNum,indTresh} = results(1).threshold.thresh;
            
        end % for {indTresh}
    end % for {indTime}
end % for {indAR}


% Compute information criterion
% _____________________________

infocrit = struct( ...
	'AIC',	NaN(nARmax,nARmax+isNum,nThreshmax), ...	% AIC
	'AICc', NaN(nARmax,nARmax+isNum,nThreshmax), ...	% AICc
	'SBIC', NaN(nARmax,nARmax+isNum,nThreshmax), ...	% SBIC
	'HQIC', NaN(nARmax,nARmax+isNum,nThreshmax));   	% HQIC
    
for indAR = 1:nARmax
    for indTime = 1-isNum:indAR
        for indTresh = 1:nThreshmax
            if strcmpi(options_tvarx.optimcrit,'loglik')
                temp = infocriterion('logL',data{indAR,indTime+isNum,indTresh},numParams(indAR,indTime+isNum,indTresh),numObs(indAR,indTime+isNum,indTresh));
            elseif strcmpi(options_tvarx.optimcrit,'logdet')
                temp = infocriterion('ResMulti',data{indAR,indTime+isNum,indTresh},numParams(indAR,indTime+isNum,indTresh),numObs(indAR,indTime+isNum,indTresh));
            end
            infocrit.AIC(indAR,indTime+isNum,indTresh)  = temp.AIC;     % AIC
            infocrit.AICc(indAR,indTime+isNum,indTresh) = temp.AICc;    % AICc
            infocrit.SBIC(indAR,indTime+isNum,indTresh) = temp.SBIC;    % SBIC
            infocrit.HQIC(indAR,indTime+isNum,indTresh) = temp.HQIC;    % HQIC
        end % for {indTresh}
    end % for {indTime}
end % for {indAR}

[~,nARAIC]  = mdmin(infocrit.AIC);
[~,nARAICc] = mdmin(infocrit.AICc);
[~,nARSBIC] = mdmin(infocrit.SBIC);
[~,nARHQIC] = mdmin(infocrit.HQIC);
ic = [nARAIC;nARAICc;nARSBIC;nARHQIC];


% Choose information criterion
% ____________________________

switch lower(infoc)
    case 'aic'              % Akaike Criterion - Specification selected
        icVar = ic(1,:);
    case 'aicc'             % Corrected Akaike Criterion - Specification selected
        icVar = ic(2,:);
    case 'sbic'             % Schwarz-Bayesian Criterion - Specification selected
        icVar = ic(3,:);         
    case 'hqic'             % Hannan-Quinn Criterion - Specification selected
        icVar = ic(4,:);         
    otherwise
        error('tvarxic:InvalidInputs','tvarxic: infoc must be AIC, AICc, SBIC of HQIC');
end % switch {infoc}

% Adjust for transition variable
icVar(2) = icVar(2) - isNum;

% Add the threshold information, if needed
if nThreshmax == 1
    icVar = [icVar 1];
end


% Plot and Print
% ______________

switch graph
    case 0
        % No graph
    case 1
        fprintf('MODEL SELECTION');
        fprintf('\nOptimal lag, time delay and threshold number for model are %d, %d and %d, according to a %s criterion \n',icVar(1),icVar(2),icVar(3),infoc);
end % switch {graph}

end % function {varxic}


% ---------------------------------
function [y,nARmax,nThreshmax,transVar] = checkInput(y,nARmax,nThreshmax,transVar)
% checkInput: Local function to check the validity of required inputs

% y
validateattributes(y,{'numeric'},{'2d','real'},'tvarxic','y',1);

% nARmax
validateattributes(nARmax,{'numeric'},{'scalar','integer','positive','finite'},'tvarxic','nARmax',2);

% nThreshmax
validateattributes(nThreshmax,{'numeric'},{'scalar','integer','positive','>=',1,'<=',2},'tvarxic','nThreshmax',3);   

% transVar
validateattributes(transVar,{'logical','numeric'},{},'tvarxic','transVar',4);
if islogical(transVar)                                                                                              
    validateattributes(transVar,{'logical'},{'numel',size(y,2)},'tvarxic','transVar',4);
    if sum(double(transVar)) > 1
        error('tvarx:InputError','tvarx: currently, for transVar, only a single value is implemented')
    elseif sum(double(transVar)) == 0
        error('tvarx:InputError','tvarx: at least one transition variable must be declared')
    end
elseif isnumeric(transVar)
    validateattributes(transVar,{'numeric'},{'real','nrows',size(y,1),'ncols',1},'tvarxic','transVar',4);
end


end % subfunction {checkInput}


% ---------------------------------
function [infoc,graph] = checkOptions(options_ic) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Get default or user-provided options
getOptions(options_ic, ...
	'infoc',    'AIC', ...  % information criterion
    'graph',    true);      % plot information criterion      

% infoc
validateattributes(infoc,{'char'},{'nrows',1},'tvarxic','options_ic.infoc',5);
validatestring(lower(infoc),{'aic','aicc','sbic','hqic'},'tvarxic','options_ic.infoc',5);

% graph
validateattributes(graph,{'logical'},{'numel',1},'varxic','options_ic.graph',5);


end % subfunction {checkOptions}

