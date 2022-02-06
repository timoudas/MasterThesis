function [backtest,thresh] = tvarxbacktest(results,ystart,yend,test,options)
% TVARXBACKTEST: Determine the prediction and forecast accuracy of a TVARX model.
%
% Syntax:
%
%       backtest = tvarxbacktest(results,ystart,yend,test,options)
%
% Description:
%
%       TVARXBACKTEST determine the forecast accuracy of a TVARX model. To do this, we perform a
%       Monte-Carlo simulation with a specified number of sample paths for each year from the start date to the most recent prior
%       year. Given the specified number of sample paths for each year, we estimate the root mean-square error (RMSE) between
%       subsequent realizations and forecasts over the time horizon.
%
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by TVARX function
%       ystart      -   [integer]   scalar, year to start back testing (4-digit, >= beg_yr in results.calstruct)
%       yend        -   [integer]	scalar, year to end back testing (4-digit, <= beg_yr in results.calstruct)
%       test        -   [char]      string, determine the test performed
%                                       'in'    for dynamic prediction (in-sample)
%                                       'out'   for dynamic forecast (out-of-sample) 
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       keepThresh  -   [logical]   logical, keep the optimal threshold for all backtest loop   
%                                       true    keep the optimal threshold for all backtest loop (default)
%       	                            false	do not keep the optimal threshold for all backtest loop, estimate a new one for
%                                               each iteration of the backtest
%       goodfit     -   [char]      string, the goodness of fit
%                                       'MSE'       for mean squared error (mse)  
%                                       'NMSE'      for normalised mean squared error (nmse)
%                                       'RMSE'      for root mean squared error (rmse) (default)
%                                       'NRMSE'     for normalised root mean squared error (nrmse)
%                                       'MSRE'      for mean squared relative error (msre)
%                                       'RMSRE'     for root mean squared relative error (rmsre)
%                                       'MSPE'      for mean squared percentage error (mspe)
%                                       'RMSPE'     for root mean squared percentage error (rmspe)
%                                       'MAE'       for mean absolute error (mae)
%                                       'MARE'      for mean absolute relative error (mare)
%                                       'MAPE'      for mean absolute percentage error (mape)
%                                       'CC'        for coefficient of correlation (cc)
%                                       'CD'        for coefficient of determination (cd)
%                                       'CE'        for coefficient of efficiency (ce)
%                                       'MAXAE'     for maximum absolute error (maxae)
%                                       'MAXARE'    for maximum absolute relative error (maxare)
%       numPaths	-	[integer]	scalar, number of sample paths for each monte-carlo simulation
%                                   	Default: numPaths = 500
%       numPer      -	[integer]	scalar, number of simulation periods for each monte-carlo simulation
%                                       Default: numPer = 4
%       vnames      -   [char]      string, vector of variable names (only relevant when print = true). If 
%       	                        print = true and no vnames is provided, a string when default variables names will
%       	                        generated
%                                       e.g.	vnames = char('y1','y2','x1','x2');(NOTE: don't bother with fixed width using char.) 
%
% Output Arguments:
%
%       backtest	-	[double]	((yend-ystart)*freq)-by-numVars-by-numPaths array of root mean-square error
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, gfit, tsprint, tvarx, tvarxresid, tvarxsim
%
% References:
%
%       Chan, K. S. (1993). Consistency and limiting distribution of the least squares estimator of a threshold 
%       autoregressive model. The Annals of Statistics, 21, 520-533.
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

narginchk(4,5); nargoutchk(1,2);

[results,ystart,yend,test] = checkInput(results,ystart,yend,test);


% Options
% _______

if nargin < 5
    options = struct.empty;
elseif nargin == 5
	validateattributes(options,{'struct'},{},'tvarxbacktest','options',5);
end

[keepThresh,goodfit,numPer,numPaths,vnames] = checkOptions(options,results);


% Check provided exogenous inputs, if any
% _______________________________________

% Define sizes and constants
numY      = results(1).sizes.numY;
timeDelay = results(1).transition.timeDelay;
%numExo    = results(1).sizes.numExo;
numLags   = results(1).sizes.numLags;
numThresh = results(1).sizes.numThresh;
freq      = results(1).calstruct.freq;
%isCons    = results(1).constant;
%isTrend   = results(1).trend;

% Check options provided, if any
isOptions = ~isempty(results(1).options);

% Check exogenous inputs provided, if any
isExo = ~isempty(results(1).beta.exo);
switch isExo
	case 0  % No exogenous variables
    	% Nothing
	case 1  % Exogenous variables
        exo = results(1).options.exo;
end % switch {isExo}


% Check if we keep the optimal threshold or do we estimate it for each iteration within the loop 
% ______________________________________________________________________________________________

if isfield(results(1).options,'gamm')
    if ~isempty(results(1).options.gamm)
        isGamma = true;
    else
        isGamma = false;
    end
else  
	isGamma = false;
end

switch keepThresh
    case 1  % Keep the optimal threshold
        switch isGamma
            case 1  % Optimal threshold already determined exogenously
                % Nothing
            case 0  % Optimal threshold now determined exogenously
                results(1).options.gamm = results(1).threshold.thresh;
        end
    case 0  % Estimate new thresholds for each iteration
        switch isGamma
            case 1  % Optimal threshold already determined exogenously
                error('tvarxbacktest:OptionError','tvarxbacktest: Optimal threshold already determined exogenously, you cannot set options.keepThresh to false')
            case 0  % Estimate new thresholds for each iteration
                % Nothing
        end
end % switch {checkExo}


% Backtest
% ________

backtest = NaN(yend-ystart+1,numY);
thresh = NaN(yend-ystart+1,numThresh);

iter = 1;
switch lower(test)
    case 'in'   % Dynamic prediction (in-sample)
        
        for yy = ystart:yend
	
            % Find the position of the last date
            EndDate = ical(yy,freq,results(1).calstruct);

            % Set up model and estimate coefficients
            switch isOptions
                case 0  % Empty options
                    if islogical(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar,timeDelay);
                    elseif isnumeric(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar(1:EndDate),timeDelay);
                    end
                case 1  % Options
                    switch isExo
                        case 0
                            % Nothing
                        case 1
                            results(1).options.exo = exo(1:EndDate,:);
                    end
                    if islogical(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar,timeDelay,results(1).options);
                    elseif isnumeric(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar(1:EndDate),timeDelay,results(1).options);
                    end
            end
            
            % Thresholds value
            thresh(iter,:) = btresults(1).threshold.thresh;
            
            % Do predictions
            W = tvarxresid(btresults,struct('simul','nonlinear','method','montecarlo','numPer',numPer,'numPaths',numPaths));
            switch isExo
                case 0
                    if islogical(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar,results(1).data.y(EndDate-numLags+1-numPer:EndDate-numPer,:));
                    elseif isnumeric(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar(EndDate+1-numPer:EndDate),results(1).data.y(EndDate-numLags+1-numPer:EndDate-numPer,:));
                    end
                case 1
                    if islogical(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar,results(1).data.y(EndDate-numLags+1-numPer:EndDate-numPer,:),exo(EndDate+1-numPer:EndDate,:));
                    elseif isnumeric(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar(EndDate+1-numPer:EndDate),results(1).data.y(EndDate-numLags+1-numPer:EndDate-numPer,:),exo(EndDate+1-numPer:EndDate,:));
                    end
            end
            eYsim = mean(Ysim,3);

            % Assess Prediction Quality
            for i = 1:numY
                gf = gfit(results(1).data.y(EndDate+1-numPer:EndDate,i),eYsim(1:numPer,i),upper(goodfit));
                backtest(iter,i) = gf.(upper(goodfit));
            end
            iter = iter + 1;		

        end % for {yy}
        
    case 'out'  % Dynamic forecast (out-of-sample)

        for yy = ystart:yend
	
            % Find the position of the last date
            EndDate = ical(yy,freq,results(1).calstruct);

            % Set up model and estimate coefficients
            switch isOptions
                case 0  % Empty options
                    if islogical(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar,timeDelay);
                    elseif isnumeric(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar(1:EndDate),timeDelay);
                    end
                case 1  % Options
                    switch isExo
                        case 0
                            % Nothing
                        case 1
                            results(1).options.exo = exo(1:EndDate,:);
                    end
                    if islogical(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar,timeDelay,results(1).options);
                    elseif isnumeric(results(1).transition.transVar)
                        btresults = tvarx(results(1).data.y(1:EndDate,:),numLags,numThresh,results(1).transition.transVar(1:EndDate),timeDelay,results(1).options);
                    end
            end

            % Thresholds value
            thresh(iter,:) = btresults(1).threshold.thresh;

            % Do forecasts
            W = tvarxresid(btresults,struct('simul','nonlinear','method','montecarlo','numPer',numPer,'numPaths',numPaths));
            switch isExo
                case 0
                    if islogical(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar,results(1).data.y(EndDate-numLags+1:EndDate,:));
                    elseif isnumeric(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar(EndDate+1:EndDate+numPer),results(1).data.y(EndDate-numLags+1:EndDate,:));
                    end
                case 1
                    if islogical(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar,results(1).data.y(EndDate-numLags+1:EndDate,:),exo(EndDate+1:EndDate+numPer,:));
                    elseif isnumeric(results(1).transition.transVar)
                        Ysim = tvarxnlsim(btresults,W,results(1).transition.transVar(EndDate+1:EndDate+numPer),results(1).data.y(EndDate-numLags+1:EndDate,:),exo(EndDate+1:EndDate+numPer,:));
                    end
            end
            eYsim = mean(Ysim,3);

            % Assess Forecast Quality
            for i = 1:numY
                gf = gfit(results(1).data.y(EndDate+1:EndDate+numPer,i),eYsim(1:numPer,i),upper(goodfit));
                backtest(iter,i) = gf.(upper(goodfit));
            end
            iter = iter + 1;		

        end % for {yy}

end

% Print results
% _____________

switch freq
    case 1  % Annual
        switch lower(test)
            case 'in'
                fprintf('%s of Actual vs Model Prediction with Horizon of %d Year(s)\n',upper(goodfit),numPer);
            case 'out'
                fprintf('%s of Actual vs Model Forecast with Horizon of %d Year(s)\n',upper(goodfit),numPer);
        end
        
	case 4  % Quarterly
        switch lower(test)
            case 'in'
                fprintf('%s of Actual vs Model Prediction with Horizon of %d Quarter(s)\n',upper(goodfit),numPer);
            case 'out'
                fprintf('%s of Actual vs Model Forecast with Horizon of %d Quarter(s)\n',upper(goodfit),numPer);
        end

	case 12 % Monthly
        switch lower(test)
            case 'in'
                fprintf('%s of Actual vs Model Prediction with Horizon of %d Month(s)\n',upper(goodfit),numPer);
            case 'out'
                fprintf('%s of Actual vs Model Forecast with Horizon of %d Month(s)\n',upper(goodfit),numPer);
        end
        
end % switch {results.calstruct.freq}
tsprint(backtest,cal(ystart,1,1,yend-ystart+1),struct('vnames',vnames(1:numY,:)))


end % function {tvarxbacktest}


% ---------------------------------
function [results,ystart,yend,test] = checkInput(results,ystart,yend,test)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxbacktest','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxbacktest:InputError','tvarxbacktest: The results structure provided must be returned by TVARX function')
end
if isempty(results(1).calstruct)
    error('tvarxbacktest:InputError','tvarxbacktest: To perform backtest accuracy, you must provide a structure returned by cal.m to your TVARX model')
end

% ystart
validateattributes(ystart,{'numeric'},{'real','finite','>=',results(1).calstruct.begin_yr},'tvarxbacktest','ystart',2);

% yend
validateattributes(yend,{'numeric'},{'real','finite','>=',results(1).calstruct.begin_yr},'tvarxbacktest','yend',3);

% test
validateattributes(test,{'char'},{},'tvarxbacktest','test',4);
validatestring(lower(test),{'in' 'out'},'tvarxbacktest','test',4);


end % subfunction {checkInput}


% ---------------------------------
function [keepThresh,goodfit,numPer,numPaths,vnames] = checkOptions(options,results) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes and constants
numY      = results(1).sizes.numY;
numExo    = results(1).sizes.numExo;
freq      = results(1).calstruct.freq;
isCons    = results(1).constant;
isTrend   = results(1).trend;

% Get default or user-provided options
getOptions(options, ...
	'keepThresh',   true, ...   % keep the estimated threshold for all the backtest
    'goodfit',      'RMSE', ...	% the goodness of fit
	'numPaths',     500, ...	% number of sample paths for each monte-carlo simulation
    'numPer',       freq, ...	% number of simulation periods for each monte-carlo simulation
    'vnames',       char.empty(numY+numExo-double(isCons)-double(isTrend),0));	% vector of variable names

% keepThresh
validateattributes(keepThresh,{'logical'},{'numel',1},'tvarxbacktest','keepThresh',5);

% goodfit
validateattributes(goodfit,{'char'},{},'tvarxbacktest','options.goodfit',5);
validatestring(upper(goodfit),{'MSE' 'NMSE' 'RMSE' 'NRMSE' 'MSRE' 'RMSRE' 'MSPE' 'RMSPE' 'MAE' 'MARE' 'MAPE' 'CC' 'CD' 'CE' 'MAXAE' 'MAXARE'}, ...
    'tvarxbacktest','options.goodfit',5);

% numPer
validateattributes(numPer,{'numeric'},{'real','finite','scalar','positive'},'tvarxbacktest','options.numPer',5);

% numPaths
validateattributes(numPaths,{'numeric'},{'real','finite','scalar','positive'},'tvarxbacktest','options.numPaths',5);

% vnames
validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'tvarxbacktest','options.vnames',5);


end % subfunction {checkOptions}

