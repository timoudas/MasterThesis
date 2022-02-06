function results = tsdate(begin_yr,begin_per,freq,numObs)
% TSDATE: Produce a time-series date string for an observation number
%         given beginning year, beginning period and frequency of the data.
%
% Syntax: 
%
%       results = tsdate(begin_yr,begin_per,freq,numObs), or
%       tsdate(begin_yr,begin_per,freq,numObs), or
%       tsdate(calstruct,numObs)
%
% Description:
%
%       Produce a time-series date string for an observation number given beginning year, beginning period and frequency of the data.
%
%       Examples:   tsdate(1974,1,12,13) would print: Jan75
%                   tsdate(1974,1,4,13) would print: Q1-77
%                   tsdate(1974,1,1,13) would print 1986
%                   out = tsdate(1974,1,12,13) would return a string `Jan75'
%
%                   cstr = cal(1974,1,12)
%                   tsdate(cstr,13) would print Jan75
%
% Input Arguments:
%
%       begin_yr    -   [integer]	scalar, the beginning year
%       begin_per	-   [integer]	scalar, the beginning period
%       freq        -   [integer]	scalar, the frequency
%                                       1 for annual
%                                       4 for quarterly
%                                       12 for monthly
%       numObs      -   [integer]   scalar, the observation number
%       or
%       calstruct	-   [struct]    a structure returned by cal()
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       results -   [char]  string of dates or simply displays the date associated with numObs (observation number)
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
%       (c) James P. LeSage, Dept of Economics
%           University of Toledo
%           2801 W. Bancroft St,
%           Toledo, OH 43606
%           jpl@jpl.econ.utoledo.edu


% Input and output arguments checking
narginchk(2,4); nargoutchk(1,1);

switch nargin
    case 4   % Case where a cal() structure is not used
        % Nothing

    case 2   % Case where a cal() structure is used
        % Error checking
        if ~isstruct(begin_yr)
            error('tsdate:InputError','tsdate: one of the two inputs must be a structure from cal');
        end
        numObs = begin_per;   % the second argument is an observation #
                            % the first argument is a structure returned by cal
        freq = begin_yr.freq;
        begin_per = begin_yr.begin_per;
        begin_yr = begin_yr.begin_yr;

    otherwise
        error('tsdate:InputError','tsdate: wrong number of arguments');

end % switch {nargin}

% Error check for 4-digit years on input
strsize = size(num2str(begin_yr),2);
if strsize ~= 4
    error('tsdate:InputError','tsdate: input a 4-digit year');
end % if {strsize ~= 4}

% Produce a time-series date string 
switch freq
     case 1	% Case of annual series 
        if begin_per > 1
            error('tsdate:InputError','tsdate: wrong begin_per argument');
        end
        ydigit = 'yyyy';  
        d = datenum(begin_yr,12*numObs,1);
        if nargout == 0
            fprintf('%6s \n',datestr(d,ydigit));
        else
            results = datestr(d,ydigit);
        end
  
	case 4	% Case of quarterly series
        if begin_per > 4
            error('tsdate:InputError','tsdate: wrong begin_per argument');
        end      
      
        ydigit = 'QQ-YY';  
        d = datenum(begin_yr,begin_per*3+3*numObs-3,1);
        if nargout == 0
            fprintf('%6s \n',datestr(d,ydigit));
        else
            results = datestr(d,ydigit);
        end

	case 12	% Case of monthly series
        if begin_per > 12
            error('tsdate:InputError','tsdate: wrong begin_per argument');
        end      
            
        ydigit = 'mmmyy';  
        d = datenum(begin_yr,begin_per+numObs-1,1);
        if nargout == 0
            fprintf('%6s \n',datestr(d,ydigit));
        else
            results = datestr(d,ydigit);
        end
  
	otherwise % How did we get here?
        disp('frequency unknown to tsdate');

end % switch {freq}

end % function {tsdate}

