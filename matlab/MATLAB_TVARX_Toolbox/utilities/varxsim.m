function Y = varxsim(results,W,Y0,X)
% VARXSIM: Simulates sample paths of a multidimensional VARX time series process.
%
% Syntax:
%
%       Y = varxsim(results,W,Y0,X)
%
% Description:
%
%       VARXSIM simulates sample paths of a multidimensional VARX time series process.
%
% Input Arguments:
%
%       results	-   [struct]	a results structure returned by VARX function
%       W       -	[double]    numObs-by-numY-numPaths innovations data
%
% Optional Input Arguments:
%
%       Y0      -	[double]    numObs-by-numPresampY0 or numObs-by-numPresampY0-by-numPaths presampled response data. Y0 is a matrix or a 
%                               3D array. If Y0 is empty or if numPresampY0 is less than the maximum AR lag in results, presample
%                               values are padded with zeros. If numPresampY0 is greater than the maximum AR lag, the most recent samples 
%                               from the last rows of each path of Y0 are used. If Y has multiple paths, Y0 must contain either
%                               a single path (applied to all paths in Y) or at least as many paths as in Y (extra paths are ignored).
%                               	Default: Y0 = []
%       X       -	[double]    numObs-by-numX or numObs-by-numX-by-numPaths exogenous data. X is a matrix or a 3D array. If Y has 
%                               multiple paths, X must contain either a single path (applied to all paths in Y) or at least as many 
%                               paths as in y (extra paths are ignored).
%                               	Default: X = []
%
% Output Arguments:
%
%       Y       -	[double]	numObs-by-numVars-by-numPaths array of simulated paths
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       MATLAB products: MATLAB
%       Functions: none
%       Subfuntions: checkInput, checkOptions
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
%       (c) Gabriel Bruneau, 2012-2016


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(0,1);

[results,W] = checkInput(results,W);

% Define constant
[numObs,~,numPaths] = size(W);
numLags = results.sizes.numLags;
numY    = results.sizes.numY;


% Options
% _______

if nargin < 3 || (nargin >= 3 && isempty(Y0))
	Y0 = zeros(numLags,numY,numPaths);
end

if nargin < 4
    X = [];
end

[Y0,X] = checkOptions(Y0,X,results,W);


% Initialization
% ______________

Y       = zeros(size(W));
ARlag   = 1:numLags;
nAR     = max(ARlag);

isCons  = results.constant;
isTrend = results.trend;
isExo   = (~isempty(results.beta.exo) && ~isempty(X));


% Simulate data
% _____________

for indPaths = 1:numPaths
	for indObs = 1:numObs

        % Add constant term to response
        if isCons
			Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + results.beta.cons.';
        end

        % Add linear trend term to response
        if isTrend
			Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + indObs*results.beta.trend.';
        end
        
		% Add exogenous inputs to response
        if isExo
            Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (X(indObs,:,indPaths)*results.beta.exo.');
        end        

        % Add autoregressive terms to response
        for indLags = 1:numLags
            if indObs <= ARlag(indLags)
				Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (Y0(nAR-ARlag(indLags)+indObs,:,indPaths)*results.beta.lags{indLags}.');
            else
                Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + (Y(indObs-ARlag(indLags),:,indPaths)*results.beta.lags{indLags}.');
            end
        end % for {indLags}
        
        % Add current innovation to response
        Y(indObs,:,indPaths) = Y(indObs,:,indPaths) + W(indObs,:,indPaths);

	end % for {indObs}
end % for {indPaths}


end % function {varxsim}


% ---------------------------------
function [results,W] = checkInput(results,W)
% checkInput: Local function to check the validity of required inputs

% Define constant
numY = results.sizes.numY;

% results
validateattributes(results,{'struct'},{},'varxsim','results',1);
if ~strcmpi(results.meth,'varx')
    error('varxsim:InputError','varxsim: The results structure provided must be returned by VARX function')
end

% W
validateattributes(W,{'numeric'},{'real','finite','size',[NaN numY NaN]},'varxsim','W',2);


end % subfunction {checkInput}


% ---------------------------------
function [Y0,X] = checkOptions(Y0,X,results,W)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes
[~,numY,numPaths] = size(W);
numLags = results.sizes.numLags;

% Y0
validateattributes(Y0,{'numeric'},{'real','finite','nrows',numLags,'ncols',numY},'varxsim','Y0',3);
if size(Y0,3) ~= 1 && size(Y0,3) ~= numPaths
	error('varxsim:InputError','varxsim: the third dimension size of Y0 should be 1 or equal to its W counterpart.');
end
if size(Y0,3) == 1
	Y0 = repmat(Y0,[1 1 numPaths]);
end

% X
if ~isempty(X)
	validateattributes(X,{'numeric'},{'real','finite','nrows',size(W,1),'ncols',size(results.beta.exo,2)},'varxsim','X',4);
    if size(X,3) ~= 1 && size(X,3) ~= numPaths
        error('varxsim:InputError','varxsim: the third dimension size of X should be 1 or equal to its W counterpart.');
    end
    if size(X,3) == 1
        X = repmat(X,[1 1 numPaths]);
    end
end


end % subfunction {checkOptions}

