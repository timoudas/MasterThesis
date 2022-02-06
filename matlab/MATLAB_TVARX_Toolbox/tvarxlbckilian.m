function [bcresults,biasConstantfinal,biasTrendfinal,biasExofinal,biasARfinal] = tvarxlbckilian(results,options)
% TVARXLBCKILIAN: Kilian (1998) bias-correction for linear (regime-specific) impulse response functions of a TVARX(p) process.
%
% Syntax:
%
%       [bcresults,biasConstantfinal,biasTrendfinal,biasExofinal,biasARfinal] = tvarxlbckilian(results,options)
%
% Description:
%
%       Killian (1998) bias-correction for linear (regime-specific) impulse response functions of a TVARX(p) process.
%
% Input Arguments:
%
%       results             -   [struct]	a results structure returned by TVARX function
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       boot                -   [char]      string for the method to get confidence intervals
%                                               'iid'	for bootstrapped confidence intervals using an iid bootstrap
%
% Output Arguments:
%
%       bcresults           -   [struct]    a structure with the same fields as tvarx    
%       biasConstantfinal   -   [cell]      numRegimes-by-1 cell array, with each cell containing regime-specific bias of constant terms
%       biasTrendfinal      -   [cell]      numRegimes-by-1 cell array, with each cell containing regime-specific bias of trend terms
%       biasExofinal        -   [cell]      numRegimes-by-1 cell array, with each cell containing regime-specific bias of exogenous terms
%       biasARfinal         -   [cell]      numRegimes-by-1 cell array, with each cell containing regime-specific bias of AR terms
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, varx, varxstab, tvarxresid, tvarxlsim, tvarxstab
%
% References:
%
%       none
%
% Notes:
%
%       TO DO: 1) Implement other bootstrap method for bias correction: continuous bootstrap (unit root only), circular block bootstrap and 
%                 stationary bootstrap.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(5,5);

results = checkInput(results);


% Options
% _______

if nargin < 2
    options = struct.empty;
elseif nargin == 2
	validateattributes(options,{'struct'},{},'tvarxlbckilian','options',2);
end

boot = checkOptions(options);


% Check the stability of the TVAR(p)
% _________________________________

isstable = tvarxstab(results);


% Initialization
% ______________

bcresults = results;
biasConstantfinal = cell(size(results,1),1);
biasTrendfinal = cell(size(results,1),1);
biasExofinal = cell(size(results,1),1);
biasARfinal = cell(size(results,1),1);


% Compute the bias-correction
% ___________________________

% Loop
for indRegimes = 1:size(results,1)

    % For speed
    isCons  = results(indRegimes).constant;
    isTrend = results(indRegimes).trend;
    isExo   = (results(indRegimes).sizes.numExo - double(isCons) - double(isTrend) ~= 0);
    numY    = results(indRegimes).sizes.numY;
    numExo  = results(indRegimes).sizes.numExo;
    numLags = results(indRegimes).sizes.numLags;
    numObs	= results(indRegimes).sizes.numObsNaNLag;
    numRegimes = results(indRegimes).sizes.numRegimes;
    
    switch isstable(indRegimes)
        case 0
            % No adjustment to bcresults, and
            biasConstantfinal{indRegimes} = 0;
            biasTrendfinal{indRegimes} = 0; 
            biasExofinal{indRegimes} = 0; 
            biasARfinal{indRegimes} = zeros(numY,numY,numLags);

        case 1        
            % Number of bootstrap replication
            numBoot = 1000;

            % Preallocate matrices for speed
            switch isCons
                case 1
                    bootConstant = zeros(numY,numBoot);
                case 0
                    % No constant
            end

            switch isTrend
                case 1
                    bootTrend = zeros(numY,numBoot);
                case 0
                    % No linear trend
            end

            switch isExo
                case 1
                    bootExo = zeros(numY,numExo - double(isCons) - double(isTrend),numBoot);
                case 0
                    % No exogenous regressors
            end

            bootAR = zeros(numY,numY,numLags,numBoot);

            % Generate bootstrapped residuals, initial values and simulated data
            bootW = tvarxresid(results(indRegimes),struct('simul','linear','method',boot,'numPaths',numBoot));

            indicesbootY0 = ceil(rand(numBoot,numLags)*(numObs));
            bootY0 = zeros(numLags,numY,numBoot);
            for indLags = 1:numBoot
                bootY0(:,:,indLags) = results(indRegimes).data.yreg(indicesbootY0(indLags,:)',:);
            end % for {i}

            switch isExo
                case 0
                    bootY1 = tvarxlsim(results(indRegimes),bootW,{bootY0});
                case 1
                    exotemp = results(indRegimes).options.exo(numLags+1:end,:);
                    exotemp = exotemp(results(indRegimes).transition.transVarInd == indRegimes,:);
                    bootY1 = tvarxlsim(results(indRegimes),bootW,{bootY0},{exotemp});
            end

            % Bootstrap loop
            nAccept = 0;
            iAccept = [];

            results(indRegimes).options.graph = false;  % We set them to false to avoid graph and results printing at each
            results(indRegimes).options.prtr  = false;  % iteration.
            
            % Initialize wait bar
            str = sprintf('Computing Killian Bias-Correction for Regime %d (out of %d) ...',indRegimes,numRegimes);
            h = waitbar(0,str);
            warning('off','getOptions:ImproperFields')
            
            switch isExo
            	case 0
                    resultstemp = results(indRegimes);
                case 1
                    resultstemp = results(indRegimes);
                    exotemp = resultstemp.options.exo(numLags+1:end,:);
                    exotemp = exotemp(resultstemp.transition.transVarInd == indRegimes,:);
                    exotemp = [resultstemp.options.exo(1:numLags,:);exotemp]; %#ok<AGROW>
                    resultstemp.options.exo = exotemp;
            end

            for iBoot = 1:numBoot

                % Estimation
                bootresults = varx([bootY0(:,:,iBoot);bootY1{1}(:,:,iBoot)],numLags,resultstemp.options);

                % Check the stability of the VAR(p)
                bootisStable = varxstab(bootresults);

                % Extract parameters
                switch isCons
                    case 1
                        bootConstant(:,iBoot) = bootresults.beta.cons;
                    otherwise
                        % No constant
                end

                switch isTrend
                    case 1
                        bootTrend(:,iBoot) = bootresults.beta.trend;
                    otherwise
                        % No linear trend
                end

                switch isExo
                    case 1
                        bootExo(:,:,iBoot) = bootresults.beta.exo;
                    case 0
                        % No exogenous regressors
                end

                for j = 1:numLags
                    bootAR(:,:,j,iBoot) = bootresults.beta.lags{j};
                end % j

                % Keep the stable draw. We discard unstable draw.
                switch bootisStable
                    case 1
                        % Add one to the accepted draw
                        nAccept = nAccept + 1;
                        iAccept = [iAccept iBoot]; %#ok<AGROW>
                    otherwise
                        % Do not keep this draw and there is no need to replace the warning message
                end % switch {bootisstable}
                
                % Update wait bar
                waitbar(iBoot/numBoot,h)
            end % for {iBoot}
            
            % Close wait bar
            close(h)

            % Compute bias estimates
            switch isCons
                case 1
                    biasConstant = mean(bootConstant(:,iAccept) - repmat(results(indRegimes).beta.cons,1,nAccept),2);
                case 0
                    % No constant
            end

            switch isTrend
                case 1
                    biasTrend = mean(bootTrend(:,iAccept) - repmat(results(indRegimes).beta.trend,1,nAccept),2);
                case 0
                    % No linear trend
            end

            switch isExo
                case 1
                    biasExo = zeros(numY,numExo - double(isCons) - double(isTrend));
                    for j = iAccept
                        biasExotemp = bootExo(:,:,j) - results(indRegimes).beta.exo;
                        biasExo = biasExo + biasExotemp;
                    end
                    biasExo = biasExo./nAccept;
                case 0
                    % No exogenous regressors
            end

            biasAR = zeros(numY,numY,numLags);
            for indLags = 1:numLags
                for j = iAccept
                    biasARtemp = bootAR(:,:,indLags,j) - results(indRegimes).beta.lags{indLags};
                    biasAR(:,:,indLags) = biasAR(:,:,indLags) + biasARtemp;
                end % j
            end % i
            biasAR = biasAR./nAccept;

            % Compute bias-corrected coefficients
            switch isCons
                case 1
                    bcresults(indRegimes).beta.cons = results(indRegimes).beta.cons - biasConstant;
                case 0
                    % No constant
            end

            switch isTrend
                case 1
                    bcresults(indRegimes).beta.trend = results(indRegimes).beta.trend - biasTrend;
                case 0
                    % No linear trend
            end

            switch isExo
                case 1
                    bcresults(indRegimes).beta.exo = results(indRegimes).beta.exo - biasExo;
                case 0
                    % No exogenous regressors
            end

            for indLags = 1:numLags
                bcresults(indRegimes).beta.lags{indLags} = results(indRegimes).beta.lags{indLags} - biasAR(:,:,indLags);
            end % for {indLags}

            % Check the stability of the VAR(p)
            [~,biasAReigResVar] = tvarxstab(bcresults(indRegimes));

            % Adjustment for eigenvalue higher than or equal to 1
            deltaprevious = 1;
            if biasAReigResVar >= 1
                while biasAReigResVar >= 1
                    delta = deltaprevious - 0.01;
                    switch isCons
                        case 1
                            bcresults(indRegimes).beta.cons = results(indRegimes).beta.cons - biasConstant*delta;
                            biasConstantfinal{indRegimes} = biasConstant*delta;
                        case 0
                            biasConstantfinal{indRegimes} = 0;
                    end            

                    switch isTrend
                        case 1
                            bcresults(indRegimes).beta.trend = results(indRegimes).beta.trend - biasTrend*delta;
                            biasTrendfinal{indRegimes} = biasTrend*delta;
                        case 0
                            biasTrendfinal{indRegimes} = 0;
                    end

                    switch isExo
                        case 1
                            bcresults(indRegimes).beta.exo = results(indRegimes).beta.exo - biasExo*delta;
                            biasExofinal{indRegimes} = biasExo*delta;
                        case 0
                            biasExofinal{indRegimes} = 0;
                    end            

                    biasARfinal{indRegimes} = zeros(size(biasAR));
                    for indLags = 1:numLags
                        bcresults(indRegimes).beta.lags{indLags} = results(indRegimes).beta.lags{indLags} - biasAR(:,:,indLags)*delta;
                        biasARfinal{indRegimes}(:,:,indLags) = biasAR(:,:,indLags)*delta;
                    end % for {indLags}

                    % Check the stability of the VAR(p)
                    [~,biasAReigResVar] = tvarxstab(bcresults(indRegimes));

                    deltaprevious = delta;
                end
            else
                switch isCons
                    case 1
                        biasConstantfinal{indRegimes} = biasConstant;
                    case 0
                        biasConstantfinal{indRegimes} = 0;
                end

                switch isTrend
                    case 1
                        biasTrendfinal{indRegimes} = biasTrend;
                    case 0
                        biasTrendfinal{indRegimes} = 0;
                end

                switch isExo
                    case 1
                        biasExofinal{indRegimes} = biasExo;
                    case 0
                        biasExofinal{indRegimes} = 0;
                end

                biasARfinal{indRegimes} = biasAR;
            end

    end % switch {isstable}

end % for {indRegimes}


end % function {tvarxlbckilian}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxlbckilian','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxlbckilian:InputError','tvarxlbckilian: The results structure provided must be returned by TVARX function')
end


end % subfunction {checkInput}


% ---------------------------------
function boot = checkOptions(options) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Get default or user-provided options
getOptions(options, ...
    'boot', 'iid');     % nonparametric bootstrap    

% boot
validateattributes(boot,{'char'},{},'tvarxlbckilian','options.boot',2);
validatestring(lower(boot),{'iid'},'tvarxlbckilian','options.boot',2);


end % subfunction {checkOptions}

