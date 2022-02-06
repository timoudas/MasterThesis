function varxplt(results,vnames)
% VARXPLT: Plots VARX model actual vs predicted and residuals.
%
% Syntax:
%
%       varxplt(results,vnames)
%
% Description:
%
%       Plots VARX model actual vs predicted and residuals.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by an VARX function
%
% Optional Input Arguments:
%
%       vnames      -   [char]      string, vector of variable names. 
%                                       e.g.	vnames = char('y1','y2','x1','x2')       
%
% Output Arguments:
%
%       Nothing, just plots the results.
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
%           based on plt_var.m from James P. LeSage's Econometrics Toolbox:
%               James P. LeSage, Dept of Economics
%               University of Toledo
%               2801 W. Bancroft St,
%               Toledo, OH 43606
%               jlesage@spatial-econometrics.com


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(0,0);

results = checkInput(results);

% Define sizes and constants
numY    = results.sizes.numY;
numExo  = results.sizes.numExo;
numObs  = results.sizes.numObsNaNLag;
isCons  = results.constant;
isTrend = results.trend;
isExo   = ~isempty(results.beta.exo);

% Options
% _______

if nargin < 2
    vnames = char.empty(numY+numExo-double(isCons)-double(isTrend),0);
end

vnames = checkOptions(vnames,results);


% Plotting results
% ________________

switch nargin
    case 1
        for i = 1:numY
            vnames{i,1}  = char(['Variable ',num2str(i)]);
        end
    case 2
        if isExo
            validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'varxplt','vnames',2);
        else
            validateattributes(vnames,{'char'},{'nrows',numY},'varxplt','vnames',2);
        end
end

% Plots
tt = 1:numObs;
clf;
switch isempty(results.calstruct)
    case 1
        for j = 1:numY
            subplot(2,1,1)
                plot(tt,results.data.ynan(:,j),'-',tt,results.yhat(:,j),'--');
                legend('Actual','Predicted');
                title([upper(results.meth),': Actual vs. Predicted ',vnames(j,:)]);
            subplot(2,1,2)
                plot(tt,results.resid(:,j))
                legend('Residuals');
                title([upper(results.meth),': Residuals ',vnames(j,:)]);
            pause
        end
    case 0
        for j = 1:numY
            subplot(2,1,1)
                tsplot([results.data.ynan(:,j) results.yhat(:,j)],results.calstruct)
                legend('Actual','Predicted');
                title([upper(results.meth),': Actual vs. Predicted ',vnames(j,:)]);
            subplot(2,1,2)
                tsplot(results.resid(:,j),results.calstruct)
                legend('Residuals');
                title([upper(results.meth),': Residuals ',vnames(j,:)]);
            pause
        end
end % switch {isempty(results.calstruct)}


end % function {varxplt}


% ---------------------------------
function [results,vnames] = checkInput(results,vnames)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'varxplt','results',1);
if ~strcmpi(results.meth,'varx')
    error('varxplt:InputError','varxplt: The results structure provided must be returned by VARX function')
end


end % subfunction {checkInput}


% ---------------------------------
function vnames = checkOptions(vnames,results)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes and constants
numY    = results.sizes.numY;
numExo  = results.sizes.numExo;
isCons  = results.constant;
isTrend = results.trend;

% vnames
validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'varxplt','vnames',2);


end % subfunction {checkOptions}

