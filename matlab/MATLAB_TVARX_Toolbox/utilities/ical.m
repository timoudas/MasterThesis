function obs = ical(year,period,cstructure)
% ICAL: Finds observation number associated with a year and a period, given a cal() structure.
%
% Syntax: 
%
%       obs = ical(year,period,c_str)
%
% Description:
%
%       Create a time-series calendar structure variable that associates a date with an observation number.
%
%       Examples:	cstr = cal(1982,1,12) 
%                   obs = ical(1986,1,cstr) would return 48
%                   cstr = cal(1982,1,12)
%                   obs = ical(1982,1,cstr) would return 1
%
% Input Arguments:
%
%       year            -   [integer]	scalar, the year
%       period          -   [integer]	scalar, the period (<= frequency)
%       cal_struct      -   [struct]    a structure returned by cal()
%
% Optional Input Arguments:
%
%       
%
% Output Arguments:
%
%       obs             -   [integer]   scalar, the observations number
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
% ___________________________________

narginchk(3,3); nargoutchk(0,1);

% Check that user input 4-digit year
strsize = size(num2str(year),2);
if strsize ~= 4
    error('ical:InputError','ical: input a 4-digit year');
end % if {strsize ~= 4}

% Check that begin_per not negative
if period < 0
    error('ical:InputError','ical: period < 0');
end % if {period < 0}


% Find observation
% ________________

if isstruct(cstructure)
	begin_yr = cstructure.begin_yr;
	begin_per = cstructure.begin_per;
	freq = cstructure.freq;
	if period > freq % check that period not > freq
        error('ical:InputError','ical: period > freq');
	end % if {period > freq}
	if year < begin_yr % check that year not > beg_yr
        error('ical:InputError','ical: year > beg_yr');
	end % if {year < begin_yr}
	if year > begin_yr
        obs1 = length(begin_per:freq);
        obs2 = length(begin_yr+1:year-1)*freq;
        obs3 = length(1:period);
        obs = obs1+obs2+obs3;
	else
        obs = length(begin_per:period);
	end
else
    error('ical:InputError','ical: requires a structure returned by cal');
end % if {isstruct(cstructure)}


end % function {ical}

