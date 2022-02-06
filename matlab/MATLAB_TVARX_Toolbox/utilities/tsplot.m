function tsplot(y,cstruct,options)
% TSPLOT: Time-series plot with dates and labels.
%
% Syntax: 
%
%        tsplot(y,cstruct), entire series, no variable names, or
%        tsplot(y,cstruct,options)
%
% Description:
%
%       Print time-series matrix or vector with dates and column labels.
%
%       Example:    cstr = cal(1980,1,12);
%                   tsplot(y,cstr,struct('begp',13,'endp',24), would print data for 1981
%                   tsplot(y,cstr,struct('begp',ical(1981,1,cstr),'endp',ical(1981,12,cstr)), which would plot the same data for 1981
%
% Input Arguments:
%
%       y       -   [double]	matrix (or vector) of series to be printed
%       cstruct	-   [struct]	structure, returned by cal()
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       begp	-   [integer]	scalar, the beginning observation to print
%       endp 	-   [integer]	scalar, the ending period to print,
%       vnames	-   [char]      a string matrix of names for a header
%                                   ex: vnames = char('y','x1','x2','x3');
%
% Output Arguments:
%
%       Time-series plotted with dates and labels.
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       getOptions
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
%       (c) James P. LeSage, Dept of Economics
%           University of Toledo
%           2801 W. Bancroft St,
%           Toledo, OH 43606
%           jlesage@spatial-econometrics.com


% Input and Output arguments checking
% ___________________________________

narginchk(2,3); nargoutchk(0,0);

validateattributes(y,{'numeric'},{'2d','real'},'tsplot','y',1);
validateattributes(cstruct,{'struct'},{},'tsplot','cstruc',2);


% Options
% _______

if nargin < 3
    options = struct();
elseif nargin == 3
	validateattributes(options,{'struct'},{},'tsplot','options',3);
end
getOptions(options, ...
    'begp',     1, ...
    'endp',     size(y,1), ...
    'vnames',   [], ...
    'fmt',      '%10.4f');     

validateattributes(begp,{'numeric'},{'real','finite','scalar','>=',1},'tsplot','options.begp',1); 
validateattributes(endp,{'numeric'},{'real','finite','scalar','<=',size(y,1)},'tsplot','options.endp',1); 
nflag = 0;
if ~isempty(vnames) %#ok<NODEF>
    nflag = 1;
    validateattributes(vnames,{'char'},{'nrows',size(y,2)},'tsplot','options.vnames',3);
end
validateattributes(fmt,{'char'},{},'tsplot','options.fmt',3); 


% Variables names
% _______________

numVars = size(y,2);
if nflag == 0 % no variable names supplied, make some up
	vnames = [];
	for i = 1:numVars
        if i < 10
            snames = 'series  ';
            name = [snames num2str(i)];
            vnames = [vnames
                      name]; %#ok<AGROW>
        else
            snames = 'series ';
            name = [snames num2str(i)];
            vnames = [vnames
                      name]; %#ok<AGROW>
        end
	end
    vnames = vnames(2:end,:);   % remove the empty space on the first row of vnames
else
	vsize = size(vnames,1); % error check vnames argument
	if vsize ~= numVars
        error('tsplot:InputError','tsplot: wrong number of vnames to tsplot');
	end
end

fsize = 9;             % font size
numObs = size(y(begp:endp,:),1); % find nobs, nvar
    
if numObs <= 120; % provide a grid for small samples
    grid = 'on';
else
    grid = 'off';
end

freq = cstruct.freq;
  
switch freq
	case 1      % case of annual series 
        out = cal(cstruct.begin_yr,cstruct.begin_per,cstruct.freq,begp);
        begin_yr = out.begin_yr;
        yr = begin_yr:begin_yr+numObs-1;
        yrs = yr';
        ydigit = 'yyyy';  
        plot(datenum(yrs,1,1),y(begp:endp,:));
        legend(vnames);
   
	case 4      % case of quarterly series
        yrs = zeros(numObs,1);
        qtr = zeros(numObs,1);
   
        out = cal(cstruct.begin_yr,cstruct.begin_per,cstruct.freq,begp);
        begin_yr = out.begin_yr;
        %beg_qtr = out.period;
        % BUG fix suggested by Stephen Burke
        % PhD Student
        % Faculty of Commerce and Bus. Admin, Dept. of Finance
        % University of British Columbia
        if out.period == 1 
            beg_qtr = 1;
        elseif out.period == 2
            beg_qtr = 4;
        elseif out.period == 3
            beg_qtr = 7;
        else
            beg_qtr = 10;
        end
      
        for i = 1:numObs
            yrs(i,1) = begin_yr;
            qtr(i,1) = beg_qtr;
            beg_qtr = beg_qtr+3;
            if beg_qtr > 12
                begin_yr = begin_yr+1;
                beg_qtr = 1;
            end
        end
        ydigit = 'QQ-YY';  
        plot(datenum(yrs,qtr,1),y(begp:endp,:));
        legend(vnames);

	case 12     % case of monthly series
        yrs = zeros(numObs,1);
        mth = zeros(numObs,1);
        out = cal(cstruct.begin_yr,cstruct.begin_per,cstruct.freq,begp);
        begin_yr = out.begin_yr;
        beg_mth = out.begin_per;
            
        for i = 1:numObs
            yrs(i,1) = begin_yr;
            mth(i,1) = beg_mth;
            beg_mth = beg_mth+1;
            if beg_mth > 12
                begin_yr = begin_yr+1;
                beg_mth = 1;
            end
        end
        ydigit = 'mmmyy';  
        plot(datenum(yrs,mth,1),y(begp:endp,:));
        legend(vnames);

	otherwise % how did we get here?
        disp('tsplot:Frequency unknown to tsplot');

end

set(gca,'fontsize',fsize); 
set(gca,'tickdir','in');
datetick('x',ydigit);
set(gca,'GridLineStyle',':');
set(gca,'xgrid',grid);
set(gca,'xcolor','blue');
set(gca,'xcolor','k');


end % function {tsplot}

