function result = cal(begin_yr,begin_per,freq,obs)
% CAL: Create a time-series calendar structure variable that associates a date with an observation number.
%
% Syntax: 
%
%       result = cal(begin_yr,begin_per,freq,obs), or
%       result = cal(cstruc,obs)
%
% Description:
%
%       Create a time-series calendar structure variable that associates a date with an observation number.
%
% Input Arguments:
%
%       begin_yr            -   [integer]	scalar, the beginning year
%       begin_per           -   [integer]	scalar, the beginning period
%       freq                -   [integer]	scalar, the frequency
%                                           1 for annual
%                                           4 for quarterly
%                                           12 for monthly
%       or
%       cal_struct          -   [struct]    a structure returned by cal()
%
% Optional Input Arguments:
%
%       numObs              -   [integer]   scalar, the number of observations
%
% Output Arguments:
%
%       returns, a structure with the following possible fields and possible values:
%       result.begin_yr     -   [integer]	scalar, the beginning year
%       result.begin_per	-   [integer]	scalar, the beginning period
%       result.freq         -   [integer]	scalar, the frequency
%       result.obs          -   [integer]	scalar, number of observations (if input)
%       result.end_yr       -   [integer]	scalar, year for observations (if input)
%       result.end_per      -   [integer]	scalar, period for observations (if input)
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


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(1,1);

switch nargin
    case 2      % case where user has input a structure from cal()
        if ~isstruct(begin_yr)
            error('cal:InputError','cal: requires a structure as input with 2 arguments');
        else
            obs = begin_per;
            begin_per = begin_yr.begin_per;
            freq = begin_yr.freq;
            begin_yr = begin_yr.begin_yr;
        end
        % Error checking
        % Check that user input 4-digit year
        strsize = size(num2str(begin_yr),2);
        if strsize ~= 4
            error('cal:InputError','cal: input a 4-digit year to cal');
        end
        % Check that begin_per not > freq or negative
        if begin_per > freq
            error('cal:InputError','cal: begin_per > freq in cal');
        end
        if begin_per < 0
            error('cal:InputError','cal: begin_per < 0 in cal');
        end
        
    case {3,4}
        % Error checking
        % Check that user input 4-digit year
        strsize = size(num2str(begin_yr),2);
        if strsize ~= 4
            error('cal:InputError','cal: input a 4-digit year to cal');
        end
        % Check that begin_per not > freq or negative
        if begin_per > freq
            error('cal:InputError','cal: begin_per > freq in cal');
        end
        if begin_per < 0
            error('cal:InputError','cal: begin_per < 0 in cal');
        end
        
end % switch {nargin}

% Create a time-series calendar structure
switch nargin
    case 4
        result.begin_yr = begin_yr;
        result.begin_per = begin_per;
        result.freq = freq;
        result.obs = obs;
        result.end_yr = begin_yr + fix((obs+begin_per-2)/freq); % Wang Mingjian (Guanghua School of Management, Bejing University) suggested this bug fix
        tmp = rem((begin_per+obs-1),freq);
        if tmp == 0
            tmp = freq;
        end
        result.end_per = tmp;

    case 3
        result.begin_yr = begin_yr;
        result.begin_per = begin_per;
        result.freq = freq;

    case 2
        result.begin_yr = begin_yr;
        result.begin_per = begin_per;
        result.freq = freq;
        result.obs = obs;
        result.end_yr = begin_yr + fix((obs+begin_per-2)/freq);
        tmp = rem((begin_per+obs-1),freq);
        if tmp == 0
        tmp = freq;
        end
        result.end_per = tmp;

    otherwise
        error('cal:InputError','cal: wrong number of inputs');

end % switch {nargin}


end % function {cal}

