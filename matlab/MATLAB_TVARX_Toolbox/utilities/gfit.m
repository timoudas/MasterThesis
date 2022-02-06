function gf = gfit(t,p,gFitMeasure,options)
% GFIT: Computes goodness of fit for regression model.
%
% Syntax: 
%
%       gf = gfit(t,p)
%       gf = gfit(t,p,gFitMeasure) 
%       gf = gfit(t,p,gFitMeasure,options)
%
% Description:
%
%       Computes goodness of fit for regression model, from t and p, which are target and prediction vectors respectively, 
%       and returns the value for M, which is one or several error related performance metrics.
%
%       Examples:
%           gf = gfit(t,p);                         for all statistics in list returned as vector
%           gf = gfit(t,p,'RMSE');                  for root mean squared error 
%           gf = gfit(t,p,{'RMSE'});                for root mean squared error
%           gf = gfit(t,p,{'MSE' 'RMSE' 'CE'});     for mean squared error, root mean squared error, and coefficient of 
%                        |                          efficiency
%                       \|/
%           gf = [mse rmse e]
%
%           options.verbose = true;
%           gf = gfit(t,p,'ALL',options);           for all statistics in list returned as vector with information posted to 
%                                                   the command line on each statistic
%           options.verbose = true;
%           gf = gfit(t,p,{'MSE' 'RMSE' 'CE'},options);	for mean squared error, root mean squared error, and coefficient of 
%                                                       efficiency as a vector with information on each of these also posted to 
%                                                       the command line
%
% Input Arguments:
%
%       t               -	[double]	numObs-by-1 vector, matrix or vector of target values for regression model
%       p               -	[double]	numObs-by-1 vector, matrix or vector of prediction values from regression model
%       gFitMeasure     -   [char]      a string or cell array of string values representing different form of goodness of fit 
%                                       measure as follows:
%               
%               'ALL'       - calculates all the measures below
%               'MSE'       - mean squared error (mse)  
%               'NMSE'      - normalised mean squared error (nmse)
%               'RMSE'      - root mean squared error (rmse)
%               'NRMSE'     - normalised root mean squared error (nrmse)
%               'MSRE'      - mean squared relative error (msre)
%               'RMSRE'     - root mean squared relative error (rmsre)
%               'MSPE'      - mean squared percentage error (mspe)
%               'RMSPE'     - root mean squared percentage error (rmspe)
%               'MAE'       - mean absolute error (mae)
%               'MARE'      - mean absolute relative error (mare)
%               'MAPE'      - mean absolute percentage error (mape)
%               'CC'        - coefficient of correlation (cc)
%               'CD'        - coefficient of determination (cd)
%               'CE'        - coefficient of efficiency (ce)
%               'MAXAE'     - maximum absolute error (maxae)
%               'MAXARE'    - maximum absolute relative error (maxare)
%
%       where
%           numObs is the sample size
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       options.verbose	-	[logical]	indicator, verbose output, posts some text output for the chosen measures to the 
%                                       command line
%                                               true	for posting 
%                                               false	for no posting (default)
%
% Output Arguments:
%
%       gf, a structure with the following fields:
%       gf.MSE      -   [double]	scalar, mean squared error (mse)  
%       gf.NMSE     -   [double]	scalar, normalised mean squared error (nmse)
%       gf.RMSE     -   [double]	scalar, root mean squared error (rmse)
%       gf.NRMSE    -   [double]	scalar, normalised root mean squared error (nrmse)
%       gf.MSRE     -   [double]	scalar, mean squared relative error (msre)
%       gf.RMSRE    -   [double]	scalar, root mean squared relative error (rmsre)
%       gf.MSPE     -   [double]	scalar, mean squared percentage error (mspe)
%       gf.RMSPE    -   [double]	scalar, root mean squared percentage error (rmspe)
%       gf.MAE      -   [double]	scalar, mean absolute error (mae)
%       gf.MARE     -   [double]	scalar, mean absolute relative error (mare)
%       gf.MAPE     -   [double]	scalar, mean absolute percentage error (mape)
%       gf.CC       -   [double]	scalar, coefficient of correlation (cc)
%       gf.CD       -   [double]	scalar, coefficient of determination (cd)
%       gf.CE       -   [double]	scalar, coefficient of efficiency (ce)
%       gf.MAXAE    -   [double]	scalar, maximum absolute error (maxae)
%       gf.MAXARE   -   [double]	scalar, maximum absolute relative error (maxare)
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       getOptions, vecc
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
%       (c) Durga Lal Shrestha, 2004-2005 (durgals@hotmail.com)
%           Revision:   1)  Date: 2005/07/03.
%                       2)  Modified from gFit 2008/11/07 by Richard Crozier, extended to report multiple
%                           statistics at once (i.e. can pass in a cell array of strings describing
%                           desired output), added checking for NaNs, removed superfluous code, added
%                           max error code and now handles matrices by reshaping them into vectors.
%                           Also made small improvement to error reporting.
%                       3)	Added verbose option.
%                       4)  Modified from gfit2 2012/06/13 by Gabriel Bruneau, to report the multiple statistics from a cell
%                           array to a structure. Added some measures of goodness of fit.


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(0,1);

t = vec(t);    % Reshape matrices into vectors (order of data is not important, but t and p must be in the same order)
p = vec(p);
    
if length(t) ~= length(p)
	error('gfit:InvalidDataSize','gfit: Invalid data size, size of t and p must be same');
end


% Options
% _______

if nargin < 4
    options = struct();
elseif nargin == 4
    if ~isstruct(options)
        error('gfit:OptionsNotStructure','gfit: options should be a structure array');
    end
end
getOptions(options, ...
	'verbose', false);	% verbose

if nargin > 2
	% Check if gFitMeasure is cell string array or just a string
	if ~iscell(gFitMeasure) && ischar(gFitMeasure)
        if strcmp(gFitMeasure,'ALL')
            % If the string 'ALL' is passed in, all the stats are required so make the appropriate cell string array
            gFitMeasure = {'MSE' 'NMSE' 'RMSE' 'NRMSE' 'MSRE' 'RMSRE' 'MSPE' 'RMSPE' 'MAE' 'MARE' 'MAPE' 'CC' 'CD' ...
                            'CE' 'MAXAE' 'MAXARE'}; % return all measures
        else
            % Otherwise convert string to cell string array of size 1
            gFitMeasure = {gFitMeasure};
        end
    else
        % If it is a cell array of strings, check its size
        if size(gFitMeasure,2) == 1
            % If there is only one element check it is not a request for all measures
            if strcmp(char(gFitMeasure),'ALL')
                % If the string 'ALL' is passed in, all the stats are required so make the appropriate cell string array
                gFitMeasure = {'MSE' 'NMSE' 'RMSE' 'NRMSE' 'MSRE' 'RMSRE' 'MSPE' 'RMSPE' 'MAE' 'MARE' 'MAPE' 'CC' 'CD' ...
                                'CE' 'MAXAE' 'MAXARE'}; % return all measures
            end
        end
	end
else
	% Return all measures if only two inputs, nothing will be posted to the command line
	gFitMeasure = {'MSE' 'NMSE' 'RMSE' 'NRMSE' 'MSRE' 'RMSRE' 'MSPE' 'RMSPE' 'MAE' 'MARE' 'MAPE' 'CC' 'CD' ...
                    'CE' 'MAXAE' 'MAXARE'}; 
end
    

% Goodness of fit
% _______________

% Remove NaNs from the arrays, avoid modifying them if there are no NaNs to prevent reallocation of memory
if sum(isnan(t) | isnan(p))
	inds = ~isnan(t) & ~isnan(p); 
	t = t(inds);  
	p = p(inds); 
end
    
% Calculate the error
e = t-p; 

% Preallocate structure array
gf = struct('MSE', [], ...
            'NMSE', [], ...
            'RMSE', [], ...
            'NRMSE', [], ...
            'MSRE', [], ...
            'RMSRE', [], ...
            'MSPE', [], ...
            'RMSPE', [], ...
            'MAE', [], ...
            'MARE', [], ...
            'MAPE', [], ...
            'CC', [], ...
            'CD', [], ...
            'CE', [], ...
            'MAXAE', [], ...
            'MAXARE', []); 

for indFitMeasure = 1:size(gFitMeasure,2)
	switch char(gFitMeasure(indFitMeasure))
        case {'MSE','mse','1'}                          % mean squared error
            gf.MSE = mean(e.^2);                        % 0 - perfect match between output and target
            if verbose
                disp(['Mean squared error (mse): ' num2str(gf.MSE)])
            end

        case {'NMSE','nmse','2'}                        % normalised mean squared error
            gf.NMSE = mean(e.^2)/var(t);                % 0 - perfect match 
            if verbose
                disp(['Normalised mean squared error (nmse): ' num2str(gf.NMSE)])
            end

        case {'RMSE','rmse','3'}                        % root mean squared error
            gf.RMSE = sqrt(mean(e.^2));                 % 0 - perfect match     
            if verbose
                disp(['Root mean squared error (rmse): ' num2str(gf.RMSE)])
            end
            
        case {'NRMSE','nrmse','4'}                      % normalised root mean squared error
            gf.NRMSE = sqrt(mean(e.^2)/var(t));         % 0 - perfect match
            if verbose
                disp(['Normalised root mean squared error (nrmse): ' num2str(gf.NRMSE)])
            end
            
        case {'MSRE','msre','5'}                        % mean squared relative error
            gf.MSRE = mean((e./t).^2);                  % 0 - perfect match 
            if verbose
                disp(['Mean squared relative error: ' num2str(gf.MSRE)])
            end
            
        case {'RMSRE','rmsre','6'}                      % root mean squared relative error
            gf.RMSRE = sqrt(mean((e./t).^2));           % 0 - perfect match 
            if verbose
                disp(['Root mean squared relative error: ' num2str(gf.RMSRE)])
            end
            
        case {'MSPE','mspe','7'}                        % mean squared percentage error
        	gf.MSPE = mean(((e./t)*100).^2);            % 0 - perfect match 
            if verbose
                disp(['Mean squared percentage error: ' num2str(gf.MSPE)])
            end         

        case {'RMSPE','rmspe','8'}                      % root mean squared percentage error
        	gf.RMSPE = sqrt(mean(((e./t)*100).^2));     % 0 - perfect match 
            if verbose
                disp(['Root mean squared percentage error: ' num2str(gf.RMSPE)])
            end
            
        case {'MAE','mae','9'}                          % mean absolute error
            gf.MAE = mean(abs(e));                      % 0 - perfect match
            if verbose
                disp(['Mean absolute error (mae): ' num2str(gf.MAE)])
            end
            
        case {'MARE','mare','10'}                       % mean absolute relative error
            gf.MARE = mean((abs(e./t)));                % 0 - perfect match
            if verbose
                disp(['Mean  absolute relative error (mare): ' num2str(gf.MARE)])
            end
            
        case {'MAPE','mape','11'}                       % mean absolute percentage error
        	gf.MAPE = mean(abs((e./t)*100));            % 0 - perfect match 
            if verbose
                disp(['Mean absolute percentage error: ' num2str(gf.MAPE)])
            end
            
        case {'CC','cc','12'}                           % coefficient of correlation
            cf = corrcoef(t,p);                         % 1 - perfect match
            gf.CC = cf(1,2);   
            if verbose
                disp(['Coefficient of correlation (r): ' num2str(gf.CC)])
            end
            
        case {'CD','cd','13'}                           % coefficient of determination (r-squared)
            cf = corrcoef(t,p);      
            gf.CD = cf(1,2);
            gf.CD = gf.CD^2;                            % 1 - perfect match
            if verbose
                disp(['Coefficient of determination (r-squared): ' num2str(gf.CD)])
            end
            
        case {'CE','ce','14'}                           % coefficient of efficiency
            gf.CE = 1-sum(e.^2)/sum((t-mean(t)).^2);	% 1 - perfect match
            if verbose
                disp(['Coefficient of efficiency (e): ' num2str(gf.CE)])
            end
            
        case {'MAXAE','maxae','15'}                     % maximum absolute error
            gf.MAXAE = max(abs(e));                     % 0 - perfect match
            if verbose
                disp(['Maximum absolute error: ' num2str(gf.MAXAE)])
            end
            
        case {'MAXARE','maxare','16'}                   %  maximum absolute relative error
            gf.MAXARE = max(abs(e./t));                 % 0 - perfect match
            if verbose
                disp(['Maximum absolute relative error: ' num2str(gf.MAXARE)])
            end
            
        otherwise
            error('Invalid goodness of fit measure in gFitMeasure(%d):\nIt must be one of the strings {MSE NMSE RMSE NRMSE MSRE RMSRE MSPE RMSPE MAE MARE MAPE CC CD CE MAXAE MAXARE}, but actually contained ''%s''',indFitMeasure,char(gFitMeasure(indFitMeasure)))

	end % switch {gFitMeasure}
    
end % for {indFitMeasure}


end % function

