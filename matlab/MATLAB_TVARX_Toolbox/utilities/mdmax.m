function [maxval,indval] = mdmax(data)
% MDMAX: Find the largest component and its position.
%
% Syntax: 
%
%       [maxval,indval] = mdmax(data)
%
% Description:
%
%       Find the largest component and its position.
%
% Input Arguments:
%
%       data	-	[double]	multidimensional input array
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       maxval	-	[double]	the largest component of the array data
%       indval	-   [double]    the indices of the maximum values in array data
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       ind2subv (MATLAB Lightspeed Toolbox (by Tom Minka))
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
%       (c) Gabriel Bruneau, 2012


% Input and output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,2);

validateattributes(data,{'numeric'},{},'mdmax','data',1)


% Maximum
% _______

% Finds the max of data and its position, when data is viewed as a 1D array
[maxval,position] = max(data(:)); 

% Transform the index in the 1D view to 4 indices, given the size of data
indval = ind2subv(size(data),position);


end % function {mdmax}

