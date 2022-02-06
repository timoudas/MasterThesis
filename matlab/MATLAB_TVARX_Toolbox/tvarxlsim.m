function Y = tvarxlsim(results,W,Y0,X)
% TVARXLSIM: Simulates linear (regime specific) sample paths of a multidimensional TVARX time series process.
%
% Syntax:
%
%       Y = tvarxlsim(results,W,Y0,X)
%
% Description:
%
%       TVARXLSIM simulates linear (regime specific) sample paths of a multidimensional TVARX time series process.
%
% Input Arguments:
%
%       results	-   [struct]	a results structure returned by TVARX function
%       W       -	[cell]      numRegimes-by-1 cell array, with each cell containing a numObs-by-numDims-numPaths innovations data, with 
%                               numObs being specific to each threshold
%
% Optional Input Arguments:
%
%       Y0      -	[cell]      numRegimes-by-1 cell array, with each cell containing a numObs-by-numPresampY0 or 
%                               numObs-by-numPresampY0-by-numPaths presampled response data. Y0 is a matrix or a 3D array. If Y0 is  
%                               empty or if numPresampY0 is less than the maximum AR lag in results, presample values are padded with zeros. 
%                               If numPresampY0 is greater than the maximum AR lag, the most recent samples from the last rows of each path 
%                               of Y0 are used. If Y has multiple paths, Y0 must contain either a single path (applied to all paths in Y)
%                               or at least as many paths as in Y (extra paths are ignored).
%                               	Default: Y0 = []
%       X       -	[cell]      numRegimes-by-1 cell array, with each cell containing a numObs-by-numX or numObs-by-numX-by-numPaths 
%                               exogenous data. X is a matrix or a 3D array. If Y has multiple paths, X must contain either a single path 
%                               (applied to all paths in Y) or at least as many paths as in y (extra paths are ignored).
%                               	Default: X = []
%
% Output Arguments:
%
%       Y       -	[cell]      numRegimes-by-1 cell array, with each cell containing a numObs-by-numVars-by-numPaths array of simulated 
%                               paths
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
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(0,1);

[results,W] = checkInput(results,W);

% Define constant
numRegimes = size(results,1);
numObs = cell(numRegimes,1);
numDims = cell(numRegimes,1);
numPaths = cell(numRegimes,1);
numLags = cell(numRegimes,1);

for indRegimes = 1:numRegimes
    [numObs{indRegimes},numDims{indRegimes},numPaths{indRegimes}] = size(W{indRegimes});
    numLags{indRegimes} = results(indRegimes).sizes.numLags;
end


% Options
% _______

if nargin < 3 || (nargin >= 3 && isempty(Y0))
	Y0 = cell(numRegimes,1);
    for indRegimes = 1:numRegimes
        Y0{indRegimes} = zeros(numLags{indRegimes},numDims{indRegimes},numPaths{indRegimes});
    end
end

if nargin < 4
    X = cell(size(results,1),1);
	for indRegimes = 1:numRegimes
        X{indRegimes} = double.empty(size(W{indRegimes},1),0);
	end
end

[Y0,X] = checkOptions(Y0,X,W,results);


% Loop over regimes
% _________________

Y = cell(size(results,1),1);

for indRegimes = 1:size(results,1)

    numLags = results(indRegimes).sizes.numLags;
    
    % Extract relevant matrices
    Wreg  = W{indRegimes};
    Y0reg = Y0{indRegimes};
    Xreg  = X{indRegimes};
    
    % Define size
    [numObs,~,numPaths] = size(Wreg);

    % Initialization
    Yreg    = zeros(size(Wreg));
    ARlag   = 1:numLags;
    nAR     = max(ARlag);

    isConstant = results(indRegimes).constant;
    isTrend    = results(indRegimes).trend;
    isX        = (~isempty(results(indRegimes).beta.exo) && ~isempty(Xreg));

    % Simulate data
    for indPaths = 1:numPaths
        for indObs = 1:numObs

            % Add constant term to response:
            if isConstant
                Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + results(indRegimes).beta.cons';
            end

            % Add linear trend term to response:
            if isTrend
                Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + indObs*results(indRegimes).beta.trend';
            end

            % Add exogenous inputs to response:
            if isX
                Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + (Xreg(indObs,:,indPaths)*results(indRegimes).beta.exo');
            end        

            % Add autoregressive terms to response:
            for indLags = 1:numLags
                if indObs <= ARlag(indLags)
                    Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + (Y0reg(nAR-ARlag(indLags)+indObs,:,indPaths)*results(indRegimes).beta.lags{indLags}');
                else
                    Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + (Yreg(indObs-ARlag(indLags),:,indPaths)*results(indRegimes).beta.lags{indLags}');
                end
            end % for {indLags}

            % Add current innovation to response:
            Yreg(indObs,:,indPaths) = Yreg(indObs,:,indPaths) + Wreg(indObs,:,indPaths);

        end % for {indObs}
    end % for {indPaths}
    
    Y{indRegimes} = Yreg;

end % for {indRegimes}


end % function {tvarxlsim}


% ---------------------------------
function [results,W] = checkInput(results,W)
% checkInput: Local function to check the validity of required inputs

% Define constant
numDims = results(1).sizes.numY;
numRegimes = size(results,1);

% results
validateattributes(results,{'struct'},{},'tvarxlsim','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxlsim:InputError','tvarxlsim: The results structure provided must be returned by TVARX function')
end

% W
validateattributes(W,{'cell'},{'numel',numRegimes},'tvarxlsim','W',2);
for indRegimes = 1:numRegimes
    validateattributes(W{indRegimes},{'numeric'},{'real','finite','size',[NaN numDims NaN]},'tvarxlsim','W',1);
end


end % subfunction {checkInput}


% ---------------------------------
function [Y0,X] = checkOptions(Y0,X,W,results)
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define constant
numRegimes = size(results,1);
numDims = cell(numRegimes,1);
numPaths = cell(numRegimes,1);
numLags = cell(numRegimes,1);

for indRegimes = 1:numRegimes
    [~,numDims{indRegimes},numPaths{indRegimes}] = size(W{indRegimes});
    numLags{indRegimes} = results(indRegimes).sizes.numLags;
end

% Y0
validateattributes(Y0,{'cell'},{'numel',numRegimes},'tvarxlsim','Y0',3);
for indRegimes = 1:numRegimes
    validateattributes(Y0{indRegimes},{'numeric'},{'real','finite','nrows',numLags{indRegimes},'ncols',numDims{indRegimes}},'tvarxlsim','Y0',3);
    if size(Y0{indRegimes},3) ~= 1 && size(Y0{indRegimes},3) ~= numPaths{indRegimes}
        error('tvarxlsim:InputError','tvarxlsim: the third dimension size of Y0 should be 1 or equal to its W counterpart.');
    end
    if size(Y0{indRegimes},3) == 1
        Y0{indRegimes} = repmat(Y0{indRegimes},[1 1 numPaths{indRegimes}]);
    end
end

% X    
validateattributes(X,{'cell'},{'numel',numRegimes},'tvarxlsim','X',4);
for indRegimes = 1:numRegimes
    if ~isempty(X{indRegimes})
        validateattributes(X{indRegimes},{'numeric'},{'real','finite','nrows',size(W{indRegimes},1),'ncols',size(results(indRegimes).beta.exo,2)}, ...
            'tvarxlsim','X',4);
        if size(X{indRegimes},3) ~= 1 && size(X{indRegimes},3) ~= numPaths{indRegimes}
            error('tvarxlsim:InputError','tvarxlsim: the third dimension size of X should be 1 or equal to its W counterpart.');
        end
        if size(X{indRegimes},3) == 1
            X{indRegimes} = repmat(X{indRegimes},[1 1 numPaths{indRegimes}]);
        end
    end
end


end % subfunction {checkOptions}

