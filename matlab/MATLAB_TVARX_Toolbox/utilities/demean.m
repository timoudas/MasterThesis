function [datad,datam] = demean(data)
% DEMEAN:  Produce a demeaned and a mean vector or matrix.
%
% Syntax: 
%
%       [datad,datam] = demean(data)
%
% Description:
%
%       Produce a demeaned and a mean vector or matrix, ignoring NaNs. 
%
% Input Arguments:
%
%       data	-	[double]    numObs-by-numVars, a matrix or vector
%
%       where
%           numObs is the sample size
%           numVars is the number of variables
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       datad	-	[double]    numObs-by-numVars, a matrix or vector demeaned
%       datam	-	[double]    numObs-by-numVars, a matrix or vector of mean
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
%       (c) Gabriel Bruneau, 2012


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,2);

% Check if data is vector or matrix
if ~ismatrix(data)
	error('demean:InputError','demean: data must be a vector or matrix')
end


% Demean
% ______

numObs = size(data,1);

% Demeaned data
m = nanmean(data);
datad = data - repmat(m,numObs,1);

% Mean
if nargout == 2
    datam = repmat(m,numObs,1);
end


end % function {demean}

