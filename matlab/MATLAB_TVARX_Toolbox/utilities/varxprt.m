function varxprt(results,vnames,fid)
% VARXPRT: Prints VARX model output.
%
% Syntax:
%
%       varxprt(result,vnames,fid)
%
% Description:
%
%       Prints VARX model output.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by VARX function
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
%       checkInput, checkOptions
%
% References:
%
%       none
%
% Notes:
%
%       You may use varxprt(results,[],fid) to print output to a file with no vnames.
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

results = checkInput(results);

% Define sizes and constants
numY     = results.sizes.numY;
numExo   = results.sizes.numExo;
numLags  = results.sizes.numLags;
numX     = results.sizes.numX;
numObs   = results.sizes.numObs;
logL     = results.logL;
isCons   = results.constant;
isTrend  = results.trend;
isExo    = ~isempty(results.beta.exo);

% Options
% _______

if nargin < 2
    vnames = char.empty(numY+numExo-double(isCons)-double(isTrend),0);
end

if nargin < 3
    fid = 1;
end

[vnames,fid] = checkOptions(vnames,fid,results);


% Print results
% _____________

fprintf(fid,'  Model  : %s\n',results.description);
fprintf(fid,'  Standard errors without DoF adjustment (maximum likelihood)\n');
fprintf(fid,'  Autoregressive order           : %6d \n',numLags);
fprintf(fid,'  Number of observations         : %6d \n',numObs);
fprintf(fid,'  Number of endogenous variables : %6d \n',numY);
fprintf(fid,'  Number of exogenous variables  : %6d \n',numExo);
fprintf(fid,'  Number of regressors           : %6d \n',numX);
fprintf(fid,'  Log-likelihood                 : %6d \n',logL);
if ~isempty(vnames)
	% Check vnames length
	validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'varxprt','vnames',2);
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
fprintf(fid,'  %14s %14s %14s %14s\n', '______________','______________','______________','______________');

if isCons
    for i = 1:numY
        if i == 1
            tempstring = sprintf('Constant(%d)',i);
        else
            tempstring = sprintf('        (%d)',i);
        end
        fprintf(fid,'  %14s ',tempstring);
        if ~isempty(results.beta.cons)
        	fprintf(fid,'%14.4f ',results.beta.cons(i));
        else
            fprintf(fid,'%14s ',' ');
        end
        if ~isempty(results.stderr.cons)
            fprintf(fid,'%14.4f ',results.stderr.cons(i));
        else
            fprintf(fid,'%14s ',' ');
        end
        if ~isempty(results.tstat.cons)
            fprintf(fid,'%14.4f ',results.tstat.cons(i));
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
        if ~isempty(results.beta.cons)
        	fprintf(fid,'%14.4f ',results.beta.trend(i));
        else
            fprintf(fid,'%14s ',' ');
        end
        if ~isempty(results.stderr.cons)
            fprintf(fid,'%14.4f ',results.stderr.trend(i));
        else
            fprintf(fid,'%14s ',' ');
        end
        if ~isempty(results.tstat.cons)
            fprintf(fid,'%14.4f ',results.tstat.trend(i));
        else
        	fprintf(fid,'%14s ',' ');
        end
        fprintf(fid,'\n');
    end
end

if isExo
	for ii = 1:(numX - numLags*numY - double(isCons) - double(isTrend))
        for jj = 1:numY   
            if jj == 1
                tempstring = sprintf('Exo(%d)(%d)',ii,jj);
            else
                tempstring = sprintf('        (%d)',jj);
            end
            fprintf(fid,'  %14s ',tempstring);
            if ~isempty(results.beta.exo)
                fprintf(fid,'%14.4f ',results.beta.exo(jj,ii));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results.stderr.exo)
                fprintf(fid,'%14.4f ',results.stderr.exo(jj,ii));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results.tstat.exo)
                fprintf(fid,'%14.4f ',results.tstat.exo(jj,ii));
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
            if ~isempty(results.beta.lags{i})
                fprintf(fid,'%14.4f ',results.beta.lags{i}(ii,jj));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results.stderr.lags{i})
                fprintf(fid,'%14.4f ',results.stderr.lags{i}(ii,jj));
            else
                fprintf(fid,'%14s ',' ');
            end
            if ~isempty(results.tstat.lags{i})
                fprintf(fid,'%14.4f ',results.tstat.lags{i}(ii,jj));
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
        if ~isempty(results.stats.sige)
            fprintf(fid,'%14.4f ',results.stats.sige(i,j));
        else
            fprintf(fid,'%14s ',' ');
        end
        fprintf(fid,'\n');
	end
end


end % function {varxprt}
  

% ---------------------------------
function [results,vnames] = checkInput(results,vnames)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'varxprt','results',1);
if ~strcmpi(results.meth,'varx')
    error('varxprt:InputError','varxprt: The results structure provided must be returned by VARX function')
end


end % subfunction {checkInput}


% ---------------------------------
function [vnames,fid] = checkOptions(vnames,fid,results)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes and constants
numY    = results.sizes.numY;
numExo  = results.sizes.numExo;
isCons  = results.constant;
isTrend = results.trend;

% vnames
validateattributes(vnames,{'char'},{'nrows',numY+numExo-double(isCons)-double(isTrend)},'varxprt','vnames',2);

% fid
validateattributes(fid,{'numeric'},{'real','finite','nonnan','integer','scalar'},'varxprt','fid',3);


end % subfunction {checkOptions}
