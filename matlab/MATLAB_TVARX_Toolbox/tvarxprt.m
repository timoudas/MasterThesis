function tvarxprt(results,vnames,fid)
% TVARXPRT: Prints TVARX model output.
%
% Syntax:
%
%       tvarxprt(result,vnames,fid)
%
% Description:
%
%       Prints TVARX model output.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by TVARX function
%
% Optional Input Arguments:
%
%       vnames      -   [char]      string, vector of variable names
%                                       e.g.	vnames = char('y1','y2','x1','x2');(NOTE: don't bother with fixed width using char.)       
%       fid         -   [function]	file-id for printing results to a file (defaults to the MATLAB command window)
%                                       e.g.    fid = fopen('var.out','w+');
%
% Output Arguments:
%
%       Prints results on the MATLAB command window (default) or in a specified file.
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
%       You may use tvarxprt(results,[],fid) to print output to a file with no vnames.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2013
%           based on prt_var.m from James P. LeSage's Econometrics Toolbox:
%               James P. LeSage, Dept of Economics
%               University of Toledo
%               2801 W. Bancroft St,
%               Toledo, OH 43606
%               jlesage@spatial-econometrics.com
%           and vgxdisp.m from the MATLAB Econometrics Toolbox


% Input and Output arguments checking
% ___________________________________

narginchk(1,3); nargoutchk(0,0);

validateattributes(results,{'struct'},{},'tvarxprt','results',1);

switch nargin
    case 1
        vnames = [];
        fid = 1;
    case 2
        validateattributes(vnames,{'char'},{},'tvarxprt','vnames',2);
        fid = 1;
end


% Loop over regimes
% _________________

for indRegimes = 1:size(results,1)

    % Set constant
    % ____________

    numY     = results(indRegimes).sizes.numY;
    numExo   = results(indRegimes).sizes.numExo;
    numLags  = results(indRegimes).sizes.numLags;
    numX     = results(indRegimes).sizes.numX;
    numObs   = results(indRegimes).sizes.numObs;
    logL     = results(indRegimes).optim.logL;
    isCons   = results(indRegimes).constant;
    isTrend  = results(indRegimes).trend;


    % Print results
    % _____________

    fprintf(fid,'  Model  : %s\n',results(indRegimes).description);
    fprintf(fid,'  Regime : %d\n',indRegimes);
    fprintf(fid,'  Standard errors without DoF adjustment (maximum likelihood)\n');
    fprintf(fid,'  Autoregressive order           : %6d \n',numLags);
    fprintf(fid,'  Number of observations         : %6d \n',numObs);
    fprintf(fid,'  Number of endogenous variables : %6d \n',numY);
    fprintf(fid,'  Number of exogenous variables  : %6d \n',numExo);
    fprintf(fid,'  Number of regressors           : %6d \n',numX);
    fprintf(fid,'  Log-likelihood                 : %6d \n',logL);
    if ~isempty(vnames)
        % Check vnames length
        validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'tvarxprt','vnames',2);
        fprintf('  Endogenous timeseries          :\n');
        for i = 1:numY    
            fprintf(fid,'      (%d) %s ',i,vnames(i,:));
            fprintf(fid,'\n');
        end
        if ~isempty(1:(numX - numLags*numY - double(isCons) - double(isTrend)))
            fprintf('  Exogenous timeseries           :\n');
            cnt = 1;
            for i = 1:(numX - numLags*numY - double(isCons) - double(isTrend))   
                fprintf(fid,'      (%d) %s ',cnt,vnames(numY+cnt,:));
                fprintf(fid,'\n');
                cnt = cnt + 1;
            end
        end
    end
    fprintf(fid,'  %14s %14s %14s %14s\n','Parameter','Value','Std. Error','t-Statistic');
    fprintf(fid,'  %14s %14s %14s %14s\n', ...
        '______________','______________','______________','______________');

    if isCons
        for i = 1:numY
            if i == 1
                tempstring = sprintf('Constant(%d)',i);
            else
                tempstring = sprintf('        (%d)',i);
            end
            fprintf(fid,'  %14s ',tempstring);
            if ~isempty(results(indRegimes).beta.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).beta.cons(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results(indRegimes).stderr.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).stderr.cons(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results(indRegimes).tstat.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).tstat.cons(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            fprintf(fid,'\n');
        end
    end

    if isTrend
        for i = 1:numY
            if i == 1
                tempstring = sprintf('Trend(%d)',i);
            else
                tempstring = sprintf('     (%d)',i);
            end
            fprintf(fid,'  %14s ',tempstring);
            if ~isempty(results(indRegimes).beta.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).beta.trend(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results(indRegimes).stderr.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).stderr.trend(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results(indRegimes).tstat.cons)
                fprintf(fid,'%14.4f ',results(indRegimes).tstat.trend(i));
            else
                fprintf(fid,'%14s ',' ');
            end
            fprintf(fid,'\n');
        end
    end

    if ~isempty(results(indRegimes).beta.exo)
        for ii = 1:(numX - numLags*numY - double(isCons) - double(isTrend))
            for jj = 1:numY   
                if jj == 1
                    tempstring = sprintf('Exo(%d)(%d)',ii,jj);
                else
                    tempstring = sprintf('        (%d)',jj);
                end
                fprintf(fid,'  %14s ',tempstring);
                if ~isempty(results(indRegimes).beta.exo)
                    fprintf(fid,'%14.4f ',results(indRegimes).beta.exo(jj,ii));
                else
                    fprintf(fid,'%14s ',' ');
                end
                if ~isempty(results(indRegimes).stderr.exo)
                    fprintf(fid,'%14.4f ',results(indRegimes).stderr.exo(jj,ii));
                else
                    fprintf(fid,'%14s ',' ');
                end
                if ~isempty(results(indRegimes).tstat.exo)
                    fprintf(fid,'%14.4f ',results(indRegimes).tstat.exo(jj,ii));
                else
                    fprintf(fid,'%14s ',' ');
                end
                fprintf(fid,'\n');
            end
        end
    end  

    ARlag = (1:numLags)';

    for i = 1:numLags
        for ii = 1:numY
            for jj = 1:numY
                if ii == 1 && jj == 1
                    tempstring = sprintf('AR(%d)(%d,%d)',ARlag(i),ii,jj);
                else
                    tempstring = sprintf('        (%d,%d)',ii,jj);
                end
                fprintf(fid,'  %14s ',tempstring);
                if ~isempty(results(indRegimes).beta.lags{i})
                    fprintf(fid,'%14.4f ',results(indRegimes).beta.lags{i}(ii,jj));
                else
                    fprintf(fid,'%14s ',' ');
                end
                if ~isempty(results(indRegimes).stderr.lags{i})
                    fprintf(fid,'%14.4f ',results(indRegimes).stderr.lags{i}(ii,jj));
                else
                    fprintf(fid,'%14s ',' ');
                end
                if ~isempty(results(indRegimes).tstat.lags{i})
                    fprintf(fid,'%14.4f ',results(indRegimes).tstat.lags{i}(ii,jj));
                else
                    fprintf(fid,'%14s ',' ');
                end
                fprintf(fid,'\n');
            end
        end
    end

    for i = 1:numY
        for j = 1:i
            tempstring = sprintf('Q(%d,%d)',i,j);
            fprintf(fid,'  %14s ',tempstring);
            if ~isempty(results(indRegimes).stats.sige)
                fprintf(fid,'%14.4f ',results(indRegimes).stats.sige(i,j));
            else
                fprintf(fid,'%14s ',' ');
            end
            fprintf(fid,'\n');
        end
    end

	fprintf(fid,'\n');

end % for {indRegimes}


end % function {tvarxprt}
  
