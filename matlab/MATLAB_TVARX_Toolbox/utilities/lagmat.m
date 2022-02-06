function ylag = lagmat(y,lags)
% LAGMAT: Create matrix of lagged time series.
%
% Syntax: 
%
%       ylag = lagmat(y), or
%       ylag = lagmat(y,lags)
%
% Description:
%
%       Create a matrix of lagged (time-shifted) series. Positive lags correspond to delays; negative lags correspond to 
%       leads.
%
%       Example:
%
%           Y = [(1:5)' (-2:2)']	% 2 series
%           lags = [-1 0 1];        % 3 lags
%           ylag = lagmat(y,lags)   % 2*3 = 6 lagged series
%
% Input Arguments:
%
%       y       -	[double]	numObs-by-numVars, time series data. Y may be a vector or a matrix. If Y is a vector, it 
%                               represents a single series. If Y is a numObs-by-numVars matrix, it represents numObs  
%                               observations of numVars series, with observations across any row assumed to occur at the 
%                               same time. The last observation of any series is assumed to be the most recent.
%
%       where
%           numObs is the number of observation
%           numVars is the number of series
%
% Optional Input Arguments:
%
%       lags	-	[integer]	vector of integer delays or leads, of length numLags, applied to each series in Y. The first 
%                               lag is applied to all series in Y, then the second lag is applied to all series in Y, and
%                               so forth. To include an unshifted copy of a series in the output, use a zero lag.
%                                   Default = 1 (the first lag only)
%
% Output Arguments:
%
%       ylag	-	[double]	numObs-by-(numVars*numLags) matrix of lagged versions of the series in Y. Columns of YLag 
%                               are, in order, all series in Y lagged by the first lag in lags, all series in Y lagged by the   
%                               second lag in lags, and so forth. Unspecified observations (presample and postsample data) 
%                               are padded with NaN values.          
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
%       (c) 1999-2010 The MathWorks, Inc. 
%           $Revision: 1.1.8.4 $   $Date: 2010/10/08 16:41:21 $


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(0,1);

y = checkInput(y);


% Options
% _______

if nargin < 2
    lags = 1;
end

lags = checkOptions(lags);


% Lagged and leaded matrix
% ________________________

% Assign default missing value
missingValue = NaN;     

% Cycle through the lags vector and shift the input time series. Positive lags are delays and negative lags are leads (series
% are flipped (reflected in time), adjusted for lags, and then flipped again). Series with zero lags are simply copied.

% Number of lags to apply to each time series
numLags = length(lags);
[numObs,numVars] = size(y);

% Preallocate
ylag = missingValue(ones(numObs,numVars*numLags));

for c = 1:numLags
	L = lags(c);
	columns = (numVars*(c-1)+1):c*numVars; % Columns to fill, this lag

	if L > 0        % Time delays
        ylag(L+1:end,columns) = y(1:end-L,:);
	elseif L < 0    % Time leads
        ylag(1:end+L,columns) = y(1-L:end,:);
    else            % No shifts
        ylag(:,columns) = y;
	end

end


end % function {lagmat}


% ---------------------------------
function y = checkInput(y)
% checkInput: Local function to check the validity of required inputs

% y
validateattributes(y,{'numeric'},{'2d','nonempty','real'},'lagmat','y',1)
if isvector(y)	% Check for a vector 
   y = y(:);    % Ensure a column vector
end


end % subfunction {checkInput}


% ---------------------------------
function lags = checkOptions(lags)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% lags
validateattributes(lags,{'numeric'},{'real','integer','vector'},'lagmat','lags',2);
lags = lags(:); % Ensure a column vector


end % subfunction {checkOptions}

