function tvarxplt(results,vnames)
% TVARXPLT: Plots TVARX model actual vs predicted values and residuals.
%
% Syntax:
%
%       tvarxplt(results,vnames)
%
% Description:
%
%       Plots TVARX model actual vs predicted values and residuals.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by an TVARX function
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
%       none
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
%       (c) Gabriel Bruneau, 2012-2013
%           based on plt_var.m from James P. LeSage's Econometrics Toolbox:
%               James P. LeSage, Dept of Economics
%               University of Toledo
%               2801 W. Bancroft St,
%               Toledo, OH 43606
%               jlesage@spatial-econometrics.com


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(0,0);

validateattributes(results,{'struct'},{},'tvarxplt','results',1);


% Plotting results
% ________________

for indRegimes = 1:size(results,1);

    % Set constant
    numY    = results(indRegimes).sizes.numY;
    numExo  = results(indRegimes).sizes.numExo;
    numObs  = results(indRegimes).sizes.numObsNaNLag;
    isCons  = results(indRegimes).constant;
    isTrend = results(indRegimes).trend;

    % Check if there is exogenous regressors
    isExo = ~isempty(results(indRegimes).beta.exo);

    switch nargin
        case 1
            for i = 1:numY
                vnames{i,1}  = char(['Variable ',num2str(i)]);
            end
        case 2
            if isExo
                validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'tvarxplt','vnames',2);
            else
                validateattributes(vnames,{'char'},{'nrows',numY},'tvarxplt','vnames',2);
            end
    end

    % Plots
    tt = 1:numObs;
    clf;
    switch isempty(results(indRegimes).calstruct)
        case 1
            for j = 1:numY
                subplot(2,1,1)
                    plot(tt,results(indRegimes).data.ynan(:,j),'-',tt,results(indRegimes).yhat(:,j),'--');
                    legend('Actual','Predicted');
                    title([upper(results(indRegimes).meth),' Regimes ',num2str(indRegimes),': Actual vs. Predicted ',vnames(j,:)]);
                subplot(2,1,2)
                    plot(tt,results(indRegimes).resid(:,j))
                    legend('Residuals');
                    title([upper(results(indRegimes).meth),' Regimes ',num2str(indRegimes),': Residuals ',vnames(j,:)]);
                pause
            end
        case 0
            for j = 1:numY
                subplot(2,1,1)
                    tsplot([results(indRegimes).data.ynan(:,j) results(indRegimes).yhat(:,j)],results(indRegimes).calstruct)
                    legend('Actual','Predicted');
                    title([upper(results(indRegimes).meth),' Regimes ',num2str(indRegimes),': Actual vs. Predicted ',vnames(j,:)]);
                subplot(2,1,2)
                    tsplot(results(indRegimes).resid(:,j),results(indRegimes).calstruct)
                    legend('Residuals');
                    title([upper(results(indRegimes).meth),' Regimes ',num2str(indRegimes),': Residuals ',vnames(j,:)]);
                pause
            end
    end % switch {isempty(results.calstruct)}

end % for {indRegimes}


end % function {tvarxplt}

