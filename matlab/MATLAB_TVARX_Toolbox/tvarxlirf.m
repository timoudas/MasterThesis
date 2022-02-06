function [IRF,bootIRFLowerCI,bootIRFUpperCI,bootIRFMedCI,bootIRFAll,impact] = tvarxlirf(results,options)
% TVARXLIRF: Computes linear (regime specific) impulse response functions for a TVARX.
%
% Syntax:
%
%       [IRF,bootIRFLowerCI,bootIRFUpperCI] = tvarxlirf(results,options)
%
% Description:
%
%       Computes linear (regime specific) impulse response functions and standard errors for a TVARX and provide optionally a simple graph 
%       to view results.
%
% Input Arguments:
%
%       results     -   [struct]        a results structure returned by TVARX function
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       sqrttype	- 	[char,double]	string or a numVars-by-numVars positive definite matrix. This input determines the 
%                                       type of covariance decomposition used.  If it is a string if must be one of the following:
%                                           'unscaled'	for unit (unscaled) shocks, covariance assumed to be an identity matrix
%                                           'scaled'    for scaled but uncorrelated shocks. Scale is based on estimated
%                                                       error standard deviations.
%                                           'chol'      for scaled and correlated shocks, Choleski decomposition. Scale is based on
%                                                       estimated error standard deviations. (default)
%                                           'spectral'	for scaled and correlated shocks, Spectral decomposition. Scale is based on
%                                                       estimated error standard deviations.
%                                           'qr'        for scaled and correlated shocks, QR decomposition. Scale is based on
%                                                       estimated error standard deviations.
%                                           'svd'       for scaled and correlated shocks, Singular value decomposition. Scale is based on
%                                                       estimated error standard deviations.
%                                           'long'      for scaled and correlated shocks, Long-run restrictions (Blanchard-Quah). Scale is based on estimated
%                                                       error standard deviations.
%                                           'girf'      for generalized impulse response of Peseran and Shin, which is equivalent to
%                                                       numVars reorderings of the covariance matrix where each variable is ordered
%                                                       first when computing the IRF to the shock to that variable.
%                                       If the input is a numVars-by-numVars positive definite matrix, it is used as the
%                                       covariance square root for computing the impulse response functions.
%       method      -   [char]          string for the method to the confidence intervals
%                                           'montecarlo'	for monte carlo confidence intervals (default)
%                                           'iid'          for bootstrapped confidence intervals using an iid bootstrap
%       cumsumirf	-   [logical]       logical, true to compute the cumulated IRF of the differentiated variables (in results.ydiff)
%                                           Default = false
%       expirf      -   [logical]       logical, true to delogarithmise the responses for variables declared as log (in results.ylog)
%                                           Default = false
%       numPaths    -   [double]        scalar, number of monte carlo or bootstrap replications (not relevant for method = 'asymp')
%                                           Default: 2000 (to match recommandation of Killian for small-sample bias correction)
%       numPer      -   [double]        scalar, number of period to be generated for the impulse response functions
%                                           Default: 20
%       p           -   [double]        scalar, quantile to get the lower and upper bound of the confidence interval. Could be any number,
%                                       but there is some for reference:
%                                           32  to get 68% confidence interval (approximatively equal to one-standard deviations)
%                                            5	to get 95% confidence interval (approximatively equal to two-standard deviations) (default)
%                                            1  to get 99% confidence interval
%       BCKilian	-   [logical]       true if the bias-correction of Killian (1998) is applied (only relevant for bootstrap, not 
%                                       monte carlo nor asymptotic)
%                                           true    for bias-correction
%                                           false   for nor bias-correction (default)
%       graph       -   [logical]       true if plot the IRF
%                                           true    for plotting the IRF (default)
%                                           false   for no plotting
%       vnames      -   [char]          string, vector of variable names (only relevant when print = true). If 
%                                       print = true and no vnames is provided, a string when default variables names will
%                                       generated
%                                               e.g.	vnames = char('y1','y2','x1','x2');
%
% Output Arguments:
%
%       IRF             -   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                       containingthe impulse responses where IRF(i,j,h) contains the impulse to Y(j) due to shock h at 
%                                       period i
%       bootIRFLowerCI	-   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                       containing the boostrap lower bound where BOOTIRFLOWERCI(i,j,h)
%                                       contains the (0+p/2)% bound for the impulse response of Y(j) due to shock h at period i
%       bootIRFUpperCI	-   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                       containing the boostrap upper bound where BOOTIRFUPPERCI(i,j,h)
%                                       contains the (100-p/2)% bound for the impulse response of Y(j) due to shock h at period i
%       bootIRFMedCI	-   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                       ontaining the boostrap median where BOOTIRFMEDCI(i,j,h)
%                                       contains the 50% for the impulse response of Y(j) due to shock h at period i
%       impact          -   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars matrix, structural 
%                                       matrix such that impact*epsilon = u
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, getOptions, varx, varxstab, varxresid, tvarxlsim, tvarxstab, tvarxlbckillian, tvarx, tvarxresid, varxsim, varxstab, 
%       varxcompmat
%
% References:
%
%       none
%
% Notes:
%
%       TO DO:  1) Implement Sign (with potentially Zero) restrictions
%               2) Implement AB model
%               3) Implement the asymptotic confidence intervals
%               4) Implement other bootstrap method for bias correction: continuous bootstrap (unit root only), circular block bootstrap and 
%                  stationary bootstrap.
%               5) Implement asymptotic confidence intervals
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(1,6);

results = checkInput(results);


% Options
% _______

if nargin < 2
    options = struct.empty;
elseif nargin == 2
	validateattributes(options,{'struct'},{},'tvarxlirf','options',2);    % check options
end

[sqrttype,method,cumsumirf,expirf,numPaths,numPer,p,BCKilian,graph,vnames] = checkOptions(options,results);


% For speed
% _________

isCons  = results(1).constant;
isTrend = results(1).trend;
isExo   = ~isempty(results(1).beta.exo);
numY    = results(1).sizes.numY;
numLags = results(1).sizes.numLags;
numRegimes = results(1).sizes.numRegimes;


% Kilian bias-correction 
% ______________________

switch BCKilian
    case 0      % No adjustment
        bcresults = results;
        biasConstantfinal = cell(size(results,1),1);
        biasTrendfinal = cell(size(results,1),1);
        biasExofinal = cell(size(results,1),1);
        biasARfinal = cell(size(results,1),1);
        for indRegimes = 1:size(results,1)
            biasConstantfinal{indRegimes} = 0;
            biasTrendfinal{indRegimes} = 0; 
            biasExofinal{indRegimes} = 0; 
            biasARfinal{indRegimes} = zeros(numY,numY,numLags);
        end
    case 1      % Adjustment
        if any(strcmp(method,{'iid'})) % Check if the input for adjustment is ok
            [bcresults,biasConstantfinal,biasTrendfinal,biasExofinal,biasARfinal] = tvarxlbckilian(results,struct('boot',method));
        else    % No adjustment
            warning('tvarxlirf:OptionError','tvarxlirf: BCKilian can be used only when a nonparametric bootstrap is used. The IRF computation continues without any bias-correction.')
            bcresults = results;
            biasConstantfinal = cell(size(results,1),1);
            biasTrendfinal = cell(size(results,1),1);
            biasExofinal = cell(size(results,1),1);
            biasARfinal = cell(size(results,1),1);
            for indRegimes = 1:size(results,1)
                biasConstantfinal{indRegimes} = 0;
                biasTrendfinal{indRegimes} = 0; 
                biasExofinal{indRegimes} = 0; 
                biasARfinal{indRegimes} = zeros(numY,numY,numLags);
            end
        end
end


% Choose identification scheme
% ____________________________

if ischar(sqrttype)
    impact = cell(size(results,1),1);
	
	for indRegimes = 1:size(results,1)

        % Impact matrix
        switch lower(sqrttype)
            case 'unscaled'             % Unit (unscaled) shocks, covariance assumed to be an identity matrix
                impact{indRegimes} = eye(numY);                 
            case 'scaled'               % Scaled but uncorrelated shocks.
                impact{indRegimes} = diag(sqrt(diag(bcresults(indRegimes).stats.sige)));    
            case 'chol'                 % Choleski decomposition.
                impact{indRegimes} = chol(bcresults(indRegimes).stats.sige,'lower');
            case 'spectral'             % Spectral decomposition.
                impact{indRegimes} = bcresults(indRegimes).stats.sige^(0.5);
            case 'qr'                   % QR decomposition.
                C = eye(numY);
                for indLags = 1:numLags
                    C = C + bcresults(indRegimes).beta.lags{indLags};
                end
                B0 = chol(bcresults(indRegimes).stats.sige,'lower');
                if rank(C) == numY
                    Q = qr((C\B0).');
                else
                    Q = qr((pinv(C)*B0).');
                end
                impact{indRegimes} = B0*Q;
            case 'svd'                  % Singular value decomposition.
                [V1,S,V2] = svd(bcresults(indRegimes).stats.sige);
                V = (V1+V2)/2;
                Z = diag(sqrt(diag(S)));
                impact{indRegimes} = V*Z;
            case 'girf'                 % Generalized impulse response of Pesaran and Shin
                impact{indRegimes} = zeros(numY);   
                for i = 1:numY
                    order = [i setdiff(1:numY,i)];
                    QTemp = bcresults(indRegimes).stats.sige(order,order);
                    sTemp = chol(QTemp,'lower');
                    [~,reOrder] = sort(order);
                    sTemp = sTemp(reOrder,reOrder);
                    impact{indRegimes}(:,i) = sTemp(:,i);
                end % i
            case 'long'                 % Long-run restrictions
                CompMat  = tvarxcompmat(bcresults);
                Finf_big = inv(eye(length(CompMat{indRegimes}))-CompMat{indRegimes});
                Finf     = Finf_big(1:numY,1:numY);
                D        = chol(Finf*bcresults(indRegimes).stats.sige*Finf','lower');
                impact{indRegimes} = Finf\D;
            case 'sign'                 % Sign (with potentially Zero) restrictions
                % Not yet implemented
            case 'ab'                   % Estimation of an AB model
                % Not yet implemented
        end
	end
elseif isnumeric(sqrttype)
    impact = cell(size(results,1),1);
	for indRegimes = 1:size(results,1)
        impact{indRegimes} = sqrttype;
	end    
end

% Get unscaled, scaled and orthogonalized IRF
% ___________________________________________

% Generate the shock processes
W = tvarxresid(bcresults,struct('simul','linear','method','irf','numPer',numPer));
W0 = cell(size(results,1),1);
for indRegimes = 1:size(results,1)
    W0{indRegimes} = zeros(size(W{indRegimes}));
end  

% Generate the processes with and without the shock (unscaled IRF)
Y1irf = tvarxlsim(bcresults,W);
Y0irf = tvarxlsim(bcresults,W0);

% Scaled IRF
IRF1 = cell(size(results,1),1);
IRF0 = cell(size(results,1),1);
for indRegimes = 1:size(results,1)
    IRF1{indRegimes} = zeros(numPer,numY,numY);
    IRF0{indRegimes} = zeros(numPer,numY,numY);
end

for indRegimes = 1:size(results,1)
    for i = 1:numPer
        for j = 1:numY
            ej = zeros(numY,1);
            ej(j) = 1;
            IRFtemp = reshape(Y1irf{indRegimes}(i,:,:),numY,numY) * impact{indRegimes} * ej;
            IRF1{indRegimes}(i,:,j) = reshape(IRFtemp,1,1,numY);
        end % for {j}
    end % for {i}
end % for {indRegimes}

for indRegimes = 1:size(results,1)
    for i = 1:numPer
        for j = 1:numY
            ej = zeros(numY,1);
            ej(j) = 1;
            IRFtemp = reshape(Y0irf{indRegimes}(i,:,:),numY,numY) * impact{indRegimes} * ej;
            IRF0{indRegimes}(i,:,j) = reshape(IRFtemp,1,1,numY);
        end % for {j}
    end % for {i}
end % for {indRegimes}

% Cumulated IRF
switch cumsumirf
    case 1
        for indRegimes = 1:size(results,1)
            IRF1{indRegimes}(:,results.ydiff,:) = cumsum(IRF1{indRegimes}(:,results.ydiff,:));
            IRF0{indRegimes}(:,results.ydiff,:) = cumsum(IRF0{indRegimes}(:,results.ydiff,:));
        end % for {indRegimes}
    case 0
        % Nothing
end

% Exponential IRF
switch expirf
    case 1
        for indRegimes = 1:size(results,1)
            IRF1{indRegimes}(:,results.ylog,:) = exp(IRF1{indRegimes}(:,results.ylog,:));
            IRF0{indRegimes}(:,results.ylog,:) = exp(IRF0{indRegimes}(:,results.ylog,:));
        end % for {indRegimes}
    case 0
        % Nothing
end

IRF = cellfun(@minus,IRF1,IRF0,'UniformOutput',false);	% We need to substract W0 (the no shock case) 
                                                      	% to remove the effect of the constant and 
                                                        % linear trend, if any

if nargout > 1

	bootIRFAll     = cell(size(results,1),1);
	bootIRFLowerCI = cell(size(results,1),1);
    bootIRFUpperCI = cell(size(results,1),1);
    bootIRFMedCI = cell(size(results,1),1);
	nAccept = cell(size(results,1),1);
    iAccept = cell(size(results,1),1);
    
    % Do bootstrap or monte carlo IRF to get confidence interval
    % Loop
    for indRegimes = 1:size(results,1)

        % Generate resampled or simulated residuals, initial values and simulated data
        bootW = tvarxresid(bcresults(indRegimes),struct('simul','linear','method',method,'numPaths',numPaths));   % resampled or simulated residuals

        indicesbootY0 = ceil(rand(numPaths,numLags)*(results(indRegimes).sizes.numObsNaNLag));  % initial values
        bootY0 = zeros(numLags,numY,numPaths);
        for indLags = 1:numPaths
            bootY0(:,:,indLags) = bcresults(indRegimes).data.y(indicesbootY0(indLags,:)',:);
        end % for {indLags}

        switch isExo % simulated data
            case 0
                bootY1 = tvarxlsim(bcresults(indRegimes),bootW,{bootY0});
            case 1
                exotemp = bcresults(indRegimes).options.exo(numLags+1:end,:);
                exotemp = exotemp(bcresults(indRegimes).transition.transVarInd == indRegimes,:);
                bootY1 = tvarxlsim(bcresults(indRegimes),bootW,{bootY0},{exotemp});
        end % switch {isExo}

        % Replications loop
        % Initialize IRF matrix for speed
        bootIRF = zeros(numPer,numY,numY,numPaths);
        nAccept{indRegimes} = 0;
        iAccept{indRegimes} = [];

        bcresults(indRegimes).options.graph = false;  % We set them to false to avoid graph and results printing at each
        bcresults(indRegimes).options.prtr  = false;  % iteration.

        % Initialize wait bar
        % Initialize wait bar
        str = sprintf('Computing IRF Confidence Intervals for Regime %d (out of %d) ...',indRegimes,numRegimes);
        h = waitbar(0,str);
        warning('off','getOptions:ImproperFields')

        switch isExo
            case 0
                bcresultstemp = bcresults(indRegimes);
            case 1
                bcresultstemp = bcresults(indRegimes);
                exotemp = bcresultstemp.options.exo(numLags+1:end,:);
                exotemp = exotemp(bcresultstemp.transition.transVarInd == indRegimes,:);
                exotemp = [bcresultstemp.options.exo(1:numLags,:);exotemp]; %#ok<AGROW>
                bcresultstemp.options.exo = exotemp;
        end
        
        for iBoot = 1:numPaths

            % Estimation
            bootresults = varx([bootY0(:,:,iBoot);bootY1{1}(:,:,iBoot)],numLags,bcresultstemp.options);
            
            % Apply the bias correction on coefficients
            switch isCons
                case 1
                    bootresults.beta.cons = bootresults.beta.cons - biasConstantfinal{indRegimes};
                case 0
                    % No constant
            end

            switch isTrend
                case 1
                    bootresults.beta.trend = bootresults.beta.trend - biasTrendfinal{indRegimes};
                case 0
                    % No linear trend
            end

            switch isExo
                case 1
                    bootresults.beta.exo = bootresults.beta.exo - biasExofinal{indRegimes};
                case 0
                    % No exogenous regressors
            end

            for indLags = 1:numLags
                bootresults.beta.lags{indLags} = bootresults.beta.lags{indLags} - biasARfinal{indRegimes}(:,:,indLags);
            end % for {indLags}

            % Get unscaled, scaled and orthogonalized IRF
            % Generate the processes with and without the shock (unscaled IRF)
            bootY1irf = varxsim(bootresults,W{indRegimes});
            bootY0irf = varxsim(bootresults,W0{indRegimes});

            % Scaled IRF
            bootIRF1 = zeros(numPer,numY,numY);
            for j = 1:numPer
                for k = 1:numY
                    bootek = zeros(numY,1);
                    bootek(k) = 1;
                    bootIRFtemp = reshape(bootY1irf(j,:,:),numY,numY) * impact{indRegimes} * bootek;
                    bootIRF1(j,:,k) = reshape(bootIRFtemp,1,1,numY);
                end % for {k}
            end % for {j}        

            bootIRF0 = zeros(numPer,numY,numY);
            for j = 1:numPer
                for k = 1:numY
                    bootek = zeros(numY,1);
                    bootek(k) = 1;
                    bootIRFtemp = reshape(bootY0irf(j,:,:),numY,numY) * impact{indRegimes} * bootek;
                    bootIRF0(j,:,k) = reshape(bootIRFtemp,1,1,numY);
                end % for {k}
            end % for {j}

            % Cumulated IRF
            switch cumsumirf
                case 1
                    bootIRF1(:,results.ydiff,:) = cumsum(bootIRF1(:,results.ydiff,:));
                    bootIRF0(:,results.ydiff,:) = cumsum(bootIRF0(:,results.ydiff,:));
                case 0
                    % Nothing
            end

            % Exponential IRF
            switch expirf
                case 1
                    bootIRF1(:,results.ylog,:) = exp(bootIRF1(:,results.ylog,:));
                    bootIRF0(:,results.ylog,:) = exp(bootIRF0(:,results.ylog,:));
                case 0
                    % Nothing
            end        

            bootIRF(:,:,:,iBoot) = bootIRF1 - bootIRF0;
  

            % Check the stability of the VAR(p)
            bootisstable = varxstab(bootresults);

            % Keep the stable draw. We discard unstable draw.
            switch bootisstable
                case 1
                    % Add one to the accepted draw
                    nAccept{indRegimes} = nAccept{indRegimes} + 1;
                    iAccept{indRegimes} = [iAccept{indRegimes} iBoot];
                otherwise
                    % Do not keep this draw
            end % switch {bootisstable}

            % Update wait bar
            waitbar(iBoot/numPaths,h)

        end % for {iBoot}

        % Close wait bar
        close(h)
        
        % Extract the confidence interval
        bootIRFAll{indRegimes}     = bootIRF;
        bootIRFLowerCI{indRegimes} = quantile(bootIRF(:,:,:,iAccept{indRegimes}),((0+p/2)/100),4);
        bootIRFUpperCI{indRegimes} = quantile(bootIRF(:,:,:,iAccept{indRegimes}),((100-p/2)/100),4);
        bootIRFMedCI{indRegimes}   = quantile(bootIRF(:,:,:,iAccept{indRegimes}),0.5,4);

    end % for {indRegimes}

end % if {nargout > 1}


% Plot IRF
% ________

% Simple graph to view results
switch graph
    case 1
        for indRegimes = 1:size(results,1)
            LB = zeros(numY);
            UB = zeros(numY);
            hfig = figure;
            set(hfig,'Position',[100 100 800 600])
            clf;
            for i = 1:numY
                for j = 1:numY
                    subplot(numY,numY,(i-1)*numY+j);
                    if nargout > 1
                        h = plot((1:numPer),IRF{indRegimes}(:,i,j),(1:numPer),bootIRFLowerCI{indRegimes}(:,i,j),(1:numPer),bootIRFUpperCI{indRegimes}(:,i,j),[1 numPer],[0 0]);
                        set(h(1),'LineWidth',2,'Color',[0.5  0.5  1])
                        set(h(2),'LineWidth',2,'Color',[0.25 0.25 0.5],'LineStyle',':')
                        set(h(3),'LineWidth',2,'Color',[0.25 0.25 0.5],'LineStyle',':')
                        set(h(4),'LineWidth',1,'Color',[0    0    0],'LineStyle','-')
                    else
                        h = plot((1:numPer),IRF{indRegimes}(:,i,j),[1 numPer],[0 0]);
                        set(h(1),'LineWidth',2,'Color',[0.5 0.5 1])
                        set(h(2),'LineWidth',1,'Color',[0   0   0],'LineStyle','-')
                    end
                    switch ~isempty(vnames)
                        case 0  % No names are provided
                            if i == 1
                                title(['e_' num2str(j)])
                            end
                            if j == 1
                                ylabel(['y_' num2str(i)])
                            end

                        case 1 % Names are provided
                            if i == 1
                                title(['\bf ' vnames(j,:)]);
                            end
                            if j == 1
                                ylabel(['\bf ' vnames(i,:)]);
                            end
                    end
                    axis tight
                    AX = axis;
                    spread = AX(4) - AX(3);
                    AX(4) = AX(4) + 0.05*(spread);
                    AX(3) = AX(3) - 0.05*(spread);
                    LB(i,j) = AX(3);
                    UB(i,j) = AX(4);
                end % j
            end % i
            UB = max(UB,[],2);
            LB = min(LB,[],2);
            for i = 1:numY
                for j = 1:numY
                    subplot(numY,numY,(i-1)*numY+j);
                    AX = axis;
                    AX(3) = LB(i);
                    AX(4) = UB(i);
                    axis(AX);
                end % for {j}
            end % for {i}
        end % for {indRegimes}
    case 0
        % No graph
end


end % function {tvarxlirf}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxlirf','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxlirf:InputError','tvarxlirf: The results structure provided must be returned by TVARX function')
end


end % subfunction {checkInput}


% ---------------------------------
function [sqrttype,method,cumsumirf,expirf,numPaths,numPer,p,BCKilian,graph,vnames] = checkOptions(options,results) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define constant
numY = results(1).sizes.numY;
switch ~isempty(results(1).beta.exo);
    case 1
        numExo = size(results(1).beta.exo,2);
    case 0
        numExo = 0;
end

% Get default or user-provided options
getOptions(options, ...
    'sqrttype',     'chol', ...         % type of covariance decomposition
	'method',       'montecarlo', ...	% method to get confidence intervals
    'cumsumirf',    false, ...          % cumulalted IRF
    'expirf',       false, ...          % exponential IRF
    'numPaths',     2000, ...           % number of monte carlo or bootstrap replications
    'numPer',       20, ...             % number of period to be generated for the impulse response functions
    'p',            5, ...              % quantile to get the lower and upper bound
    'BCKilian',     false, ...          % Killian bias-correction
    'graph',        true, ...           % graph the IRF
    'vnames',       char.empty(numY+numExo,0));	% vector of variable names (only relevant when graph = true)
 

% sqrttype
validateattributes(sqrttype,{'char','numeric'},{},'tvarxlirf','options.sqrttype',2);
if ischar(sqrttype)
    validatestring(lower(sqrttype),{'unscaled','scaled','chol','spectral','qr','svd','long','girf'},'tvarxlirf','options.sqrttype',2);
elseif isnumeric(sqrttype)
    validateattributes(sqrttype,{'numeric'},{'2d','real','finite','size',[numY numY]},'tvarxlirf','options.sqrttype',2);
end

% method
validateattributes(method,{'char'},{},'tvarxlirf','options.boot',2);
validatestring(lower(method),{'montecarlo' 'iid'},'tvarxlirf','options.boot',2);

% cumsumirf
validateattributes(cumsumirf,{'logical'},{'numel',1},'tvarxlirf','options.cumsumirf',2);

% expirf
validateattributes(expirf,{'logical'},{'numel',1},'tvarxlirf','options.expirf',2);

% numPaths
validateattributes(numPaths,{'numeric'},{'real','scalar','integer','finite','positive'},'tvarxlirf','options.numPaths',2);

% numPer
validateattributes(numPer,{'numeric'},{'real','scalar','integer','finite','positive'},'tvarxlirf','options.numPer',2);

% p
validateattributes(p,{'numeric'},{'real','scalar','integer','finite','>',0,'<',100},'tvarxlirf','options.p',2);

% BCKilian
validateattributes(BCKilian,{'logical'},{'numel',1},'tvarxlirf','options.BCKilian',2);

% graph
validateattributes(graph,{'logical'},{'numel',1},'tvarxlirf','options.graph',2);

% vnames
validateattributes(vnames,{'char'},{'nrows',numY+numExo},'tvarxlirf','options.vnames',2);


end % subfunction {checkOptions}


