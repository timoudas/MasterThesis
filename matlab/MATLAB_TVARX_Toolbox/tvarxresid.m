function W = tvarxresid(results,options)
% TVARXRESID: Create innovations array by simulation, bootstrap, or for impulse response function.
%
% Syntax:
%
%       resid = tvarxresid(results,options)
%
% Description:
%
%       Create innovations array by simulation, bootstrap, or for impulse response function.
%                
% Input Arguments:
%
%       results         -   [struct]    a results structure returned by TVARX function
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%           simul       -   [char]              string, indicate the type of simulation
%                                                   'linear'        for linear (regime specific) simulated path of residuals (default)
%                                                   'nonlinear'     for non linear (not regime specific) simulated path of residuals
%           method      -   [cell,char,numeric]	if a string, for the method
%                                                   'montecarlo'	for simulating innovations data (default)
%                                                   'cont'         for bootstrapping innovations data using a continuous bootstrap (unit root only) (only for simul = 'nonlinear') 
%                                                   'block'        for bootstrapping innovations data using a circular block bootstrap (only for simul = 'nonlinear')
%                                                   'stat'         for bootstrapping innovations data using a stationary bootstrap (only for simul = 'nonlinear')
%                                                   'iid'          for bootstrapping innovations data using an liid bootstrap
%                                                   'irf'          for building an innovations data matrix to be used by tvarxlirf
%                                               if cell (only for simul == linear) 
%                                                   a numRegimes-by-1 cell array of numPer-by-1 vector to impose the choice of residuals to resample
%                                               if numeric (only for simul == nonlinear) 
%                                                   an array of numPer-by-1 vector to impose the choice of residuals to resample
%           lengthBlock	-	[double,char]       a scalar or a string, average block length. Note: only relevant for 'block','stat' and 'cont' in the nonlinear 
%                                               simulation case
%                                                   if scalar, average block length
%                                                   if string, 'auto' for automatic block length selection
%           numPer      -   [double,cell]       a numRegimes-by-1 cell array, number of observations to be included in the resid array.
%                                                   Default: Number of observations adjusted for lags and NaN values
%           numPaths	-   [double]            scalar, number of paths to be included in the resid array. Note: not relevant for 'irf', where numPaths
%                                               is always equal to the number of dimensions of the endogenous variables in results.
%                                                   Default: 1
%
% Output Arguments:
%
%       W               -	[cell]      numRegimes-by-1 cell array, with each cell containing a numPer-by-numY-numPaths innovations data, with 
%                                       (potentially) numPer being specific to each threshold
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, demean, bootoptblock, bootcont, bootblock, bootstationary, bootiid
%
% References:
%
%       none
%
% Notes:
%
%       TO DO:  1) Implement other bootstrap method for linear simulation in a consistent way: continuous bootstrap (unit root only), circular block bootstrap and 
%                  stationary bootstrap.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(0,1);

results = checkInput(results);


% Options
% _______

if nargin < 2
    options = struct.empty;
elseif nargin == 2
	validateattributes(options,{'struct'},{},'tvarxresid','options',2);
end

[simul,method,lengthBlock,numPer,numPaths] = checkOptions(options,results);


% Generate innovations
% ____________________

% Define constants
numRegimes = size(results,1);

switch lower(simul)
    case 'linear'   % linear simulation
        W = cell(numRegimes,1);

        for indRegimes = 1:numRegimes
            
            % Define sizes
            numDims = results(indRegimes).sizes.numY;
            numObs  = results(indRegimes).sizes.numObsNaNLag;
    
            if iscell(method)
                % Initialize residuals array
                W{indRegimes} = zeros(numPer{indRegimes},numDims,numPaths);

                % Centered the vector of residuals
                demeanW = demean(results(indRegimes).resid);

                % Generate innovation paths
                for indPaths = 1:numPaths
                    W{indRegimes}(:,:,indPaths) = demeanW(method{indRegimes}(:,indPaths),:);
                end    
    
            elseif ischar(method)    
                switch lower(method)
                    case 'montecarlo'	% Monte Carlo simulation
                        % Initialize residuals array
                        W{indRegimes} = zeros(numPer{indRegimes},numDims,numPaths);

                        % Cholesky decomposition of the residuals covariance
                        sigechol = chol(results(indRegimes).stats.sige);

                        % Generate innovation paths
                        for indPaths = 1:numPaths
                            W{indRegimes}(:,:,indPaths) = randn(numPer{indRegimes},numDims)*sigechol;
                        end

                    case 'iid'          % IID bootstrap
                        % Initialize residuals array
                        W{indRegimes} = zeros(numPer{indRegimes},numDims,numPaths);

                        % Bootstrap time index
                        indicesbootW = bootiid((1:1:numObs)',numPaths,numPer{indRegimes});

                        % Centered the vector of residuals
                        demeanW = demean(results(indRegimes).resid);

                        % Generate innovation paths
                        for indPaths = 1:numPaths
                            W{indRegimes}(:,:,indPaths) = demeanW(indicesbootW(:,indPaths),:);
                        end

                    case 'irf'         % Linear impulse response
                        % Identity matrix, followed by zeros
                        W{indRegimes} = [reshape(eye(numDims),1,numDims,numDims); zeros(numPer{indRegimes}-1,numDims,numDims)];
                
                end
                
            end
            
        end
    
    case 'nonlinear'    % non linear simulation
        % Define sizes
        numDims = results(1).sizes.numY;
    
        % Stack innovations together
        if ischar(method) && strcmpi(method,'montecarlo')
            % Nothing
        else
            bootnumPer = 0;
            for indPer = 1:numRegimes
                bootnumPer = bootnumPer + results(indPer).sizes.numObsNaNLag;
            end
                                        
            residall = NaN(bootnumPer,numDims);
            for indInn = 1:numRegimes
                residall(results(indInn).transition.transVarInd == indInn,:) = results(indInn).resid;
            end        
        end
        
        if isnumeric(method)
        	% Initialize residuals array
            W = zeros(numPer,numDims,numPaths);
            
            % Centered the vector of residuals
            demeanW = demean(residall);

            % Generate innovation paths
            for indPaths = 1:numPaths
            	W(:,:,indPaths) = demeanW(method(:,indPaths),:);
            end    
    
        elseif ischar(method)    
            switch lower(method)
                case 'montecarlo'	% Nonlinear Monte Carlo simulation
                    % Generate innovation paths
                    W = randn(numPer,numDims,numPaths);
                
                case 'cont'         % Continuous bootstrap
                    % Initialize residuals array
                    W = zeros(numPer,numDims,numPaths);

                    % Bootstrap time index
                    indicesbootW = bootcont((1:1:bootnumPer)',lengthBlock,numPaths,numPer);

                    % Centered the vector of residuals
                    demeanW = demean(residall);

                    % Generate innovation paths
                    for indPaths = 1:numPaths
                        W(:,:,indPaths) = demeanW(indicesbootW(:,indPaths),:);
                    end

                case 'block'          % Circular block bootstrap
                    % Initialize residuals array
                    W = zeros(numPer,numDims,numPaths);

                    % Select block length
                    if strcmpi(lengthBlock,'auto')
                        blocklength = mean(bootoptblock(residall),2);
                        lengthBlock = ceil(blocklength(2));
                    end
                    
                    % Bootstrap time index
                    indicesbootW = bootblock((1:1:bootnumPer)',lengthBlock,numPaths,numPer);

                    % Centered the vector of residuals
                    demeanW = demean(residall);

                    % Generate innovation paths
                    for indPaths = 1:numPaths
                        W(:,:,indPaths) = demeanW(indicesbootW(:,indPaths),:);
                    end
                    
                case 'stat'         % Stationary bootstrap
                    % Initialize residuals array
                    W = zeros(numPer,numDims,numPaths);

                    % Select block length
                    if strcmpi(lengthBlock,'auto')
                        blocklength = mean(bootoptblock(residall),2);
                        lengthBlock = ceil(blocklength(1));
                    end
                    
                    % Bootstrap time index
                    indicesbootW = bootstationary((1:1:bootnumPer)',lengthBlock,numPaths,numPer);

                    % Centered the vector of residuals
                    demeanW = demean(residall);

                    % Generate innovation paths
                    for indPaths = 1:numPaths
                        W(:,:,indPaths) = demeanW(indicesbootW(:,indPaths),:);
                    end
                    
                case 'iid'          % IID bootstrap
                    % Initialize residuals array
                    W = zeros(numPer,numDims,numPaths);

                    % Bootstrap time index
                    indicesbootW = bootiid((1:1:bootnumPer)',numPaths,numPer);

                    % Centered the vector of residuals
                    demeanW = demean(residall);

                    % Generate innovation paths
                    for indPaths = 1:numPaths
                        W(:,:,indPaths) = demeanW(indicesbootW(:,indPaths),:);
                    end
                    
                case 'irf'        % Nonlinear impulse response
                    % Identity matrix, followed by zeros
                    W = [reshape(eye(numDims),1,numDims,numDims); zeros(numPer-1,numDims,numDims)];
                
            end % switch lower(method)
                
        end % if isnumeric(method)

end % switch lower(simul)
        

end % function {tvarxresid}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxresid','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxresid:InputError','tvarxresid: The results structure provided must be returned by TVARX function')
end


end % subfunction {checkInput}


% ---------------------------------
function [simul,method,lengthBlock,numPer,numPaths] = checkOptions(options,results) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes and constants
numRegimes = size(results,1);
numObsNaNLag = NaN(numRegimes,1);
for indRegimes = 1:numRegimes
    numObsNaNLag(indRegimes) = results(indRegimes).sizes.numObsNaNLag;
end

% Default simulation length
numPer = cell(numRegimes,1);
for indRegimes = 1:numRegimes
	numPer{indRegimes} = results(indRegimes).sizes.numObsNaNLag;
end

% Get default or user-provided options
getOptions(options, ...
	'simul',        'linear', ...       % simulation type
	'method',       'montecarlo', ...   % method of innovations data
    'lengthBlock',	'auto', ...         % number of observations to be included in the resid array
    'numPer',       numPer, ...         % number of observations to be included in the resid array
    'numPaths',     1);                 % number of paths to be included in the resid array

% simul
validateattributes(simul,{'char'},{'nrows',1},'tvarxresid','options.simul',2);
validatestring(lower(simul),{'linear' 'nonlinear'},'tvarxresid','options.simul',2);

% method
validateattributes(method,{'char','cell','numeric'},{},'tvarxresid','options.method',2); %#ok<NODEF>
if ischar(method)
    validateattributes(method,{'char'},{'nrows',1},'tvarxresid','options.method',2); 
    if strcmpi(simul,'linear')
        validatestring(lower(method),{'montecarlo' 'iid' 'irf'},'tvarxresid','options.method',2);
    elseif strcmpi(simul,'nonlinear')
        validatestring(lower(method),{'montecarlo' 'cont' 'block' 'stat' 'iid' 'irf'},'tvarxresid','options.method',2);
    end
elseif isnumeric(method)
    if strcmpi(simul,'linear')
        for indRegimes = 1:numRegimes
            validateattributes(method,{'numeric'},{'2d','real','integer','finite','positive','nonzero','size',[NaN numPaths], ...
                                            '<=',numObsNaNLag(indRegimes),'>=',1},'tvarxresid','options.method',2); % method
        end
    elseif strcmpi(simul,'nonlinear')       
        validateattributes(method,{'numeric'},{'2d','real','integer','finite','positive','nonzero','size',[NaN numPaths], ...
        	'<=',results(1).sizes.numObs,'>=',1},'tvarxresid','options.method',2); % method
    end
elseif iscell(method)
	if strcmpi(simul,'nonlinear')
        error('tvarxresid:InputError','tvarxresid: Not valid method, if simul == ''nonlinear'', a numeric method array must be provided')
	end
	validateattributes(method,{'cell'},{'numel',numRegimes},'tvarxresid','options.method',2);
    for indRegimes = 1:numel(method,1)
    	validateattributes(method{indRegimes},{'numeric'},{'2d','real','integer','finite','positive','nonzero','size',[NaN numPaths], ...
                                        '<=',numObsNaNLag(indRegimes),'>=',1},'tvarxresid','options.method',2); % method
    end
end

if isnumeric(method) && strcmpi(simul,'linear')
    methodtemp = method;
    method = cell(numRegimes,1);
	for indRegimes = 1:size(results,1)
        method{indRegimes} = methodtemp;
	end
end

% lengthBlock
validateattributes(lengthBlock,{'char','numeric'},{},'tvarxresid','options.lengthBlock',2);
if isnumeric(lengthBlock)
    % Check simulation type
    if strcmpi(simul,'linear')
        error('tvarxresid:InputError','tvarxresid: No lengthBlock required with the specified simulation type.')
    end
    % Check block length consistency
    templengthBlock = NaN(numRegimes,1);
    for indRegimes = 1:numRegimes
        templengthBlock(indRegimes) = results(indRegimes).sizes.numObsNaNLag;
    end
    minlengthBlock = min(templengthBlock);
    validateattributes(lengthBlock,{'numeric'},{'real','finite','positive','scalar','>=',1,'<=',minlengthBlock}, ...
                                        'tvarxresid','options.lengthBlock',2); % lengthBlock
elseif ischar(lengthBlock)
	validatestring(lower(lengthBlock),{'auto'},'tvarxresid','options.lengthBlock',2);
    if strcmpi(method,{'cont'}) && strcmpi(simul,'nonlinear')
        error('tvarxresid:InputError','tvarxresid: Not valid lengthBlock, if simul == ''nonlinear'' and method == ''cont'', lengthBlock must be numeric')
    end
end

% numPer
if iscell(method)
    numPer = cell(size(results,1),1);
    for indRegimes = 1:size(results,1)
        numPer{indRegimes} = size(method{indRegimes},1);
    end
end

validateattributes(numPer,{'cell','numeric'},{},'tvarxresid','options.numPer',2);
if isnumeric(numPer)
    validateattributes(numPer,{'numeric'},{'real','scalar','integer','finite','positive'},'tvarxresid','options.numPer',2);
elseif iscell(numPer)
    if strcmpi(simul,'nonlinear')
        error('tvarxresid:InputError','tvarxresid: If simul == ''nonlinear'', numPer cannot be regime specific')
    end
    validateattributes(numPer,{'cell'},{'numel',numRegimes},'tvarxresid','options.numPer',2);
    for indRegimes = 1:numRegimes
        validateattributes(numPer{indRegimes},{'numeric'},{'real','scalar','integer','finite','positive'},'tvarxresid','options.numPer',2);
    end        
end

if isnumeric(numPer) && strcmpi(simul,'linear')
    numPertemp = numPer;
    numPer = cell(numRegimes,1);
	for indRegimes = 1:size(results,1)
        numPer{indRegimes} = numPertemp;
	end
end

% numPaths
validateattributes(numPaths,{'numeric'},{'scalar','real','finite','positive'},'tvarxresid','options.numPaths',2);


end % subfunction {checkOptions}

