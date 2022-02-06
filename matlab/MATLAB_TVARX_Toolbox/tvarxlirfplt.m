function tvarxlirfplt(results,IRF,options)
% TVARXLIRFPLT: Plot the IRFs computed with TVARXLIRF.
%
% Syntax:
%
%       tvarxlirfplt(results,IRF,options)
%
% Description:
%
%       Plot the IRFs computed with TVARXLIRF.
%
% Input Arguments:
%
%       results             -   [struct]	a results structure returned by TVARX function
%       IRF                 -   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                           containingthe impulse responses where IRF(i,j,h) contains the impulse to Y(j) due to shock h at 
%                                           period i
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%           bootIRFLowerCI	-   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                           containing the boostrap lower bound where BOOTIRFLOWERCI(i,j,h)
%                                           contains the (0+p/2)% bound for the impulse response of Y(j) due to shock h at period i
%           bootIRFUpperCI	-   [cell]      numRegimes-by-1 cell array, with each cell containing a numVars-by-numVars-by-numPer matrix 
%                                           containing the boostrap upper bound where BOOTIRFUPPERCI(i,j,h)
%                                           contains the (100-p/2)% bound for the impulse response of Y(j) due to shock h at period i
%           vnames          -   [char]      string, vector of variable names
%                                               e.g. vnames = char('y1','y2','x1','x2');
%           variables       -   [logical]	a logical vector containing the variables to be plotted
%                                           Default: variables = true(1,numY)
%           shocks          -   [logical]	a logical vector containing the shocks to be plotted
%                                           Default: shocks = true(1,numY)
%           filename        -   [char]      string, name for file saving
%           fontsize        -   [double]    integer, fontsize for graph
%           suptit          -   [logical]   if true, a sup title is added above the subplot
%                                               Default: true
%
% Output Arguments:
%
%       Nothing, just plots the results.
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, checkOptions, subtightplot, maxfig, FigFont, SupTitle, export_fig
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
%           based on VARirfplot.m by Ambrogio Cesa Bianchi, May 2012 (ambrogio.cesabianchi@gmail.com)


% Input and Output arguments checking
% ___________________________________

narginchk(2,10); nargoutchk(0,0);

[results,IRF] = checkInput(results,IRF);

% Define constant
numRegimes = size(results,1);


% Options
% _______

if nargin < 3
    options = struct.empty;
elseif nargin == 3
	validateattributes(options,{'struct'},{},'tvarxirfplt','options',3);
end

[bootIRFLowerCI,bootIRFUpperCI,vnames,variables,shocks,filename,fontsize,suptit] = checkOptions(options,results,IRF);


% Loop over regimes
% _________________

for indRegimes = 1:numRegimes

    % Check optional inputs and Define some parameters
    % ________________________________________________

    % Define sizes
    [numPer,numY,numShocks] = size(IRF{indRegimes});

    % Define vector of variables
    varvec = 1:1:numY;
    varvec = varvec(variables);

    % Define vector of shocks
    shockvec = 1:1:numShocks;
    shockvec = shockvec(shocks);

    % Define the rows and columns for the subplots
    row = round(sqrt(length(varvec)));
    col = ceil(sqrt(length(varvec)));
    colvec = 1;
    for i = 1:numY 
        colvec = [colvec col*i+1]; %#ok<AGROW>
    end

    % Define a timeline
    steps = 1:1:numPer;
    x_axis = zeros(1,numPer);


    % Plot
    % ____

    subplot = @(m,n,p) subtightplot (m,n,p,[0.09 0.09]);%,[0.1 0.01],[0.1 0.01]

    figure();
    maxfig();
    % Old line
    %FigSize(22.86,16.51) % Fill a whole Letter paper size (Height:27.94cm, Width: 21.59cm), adjusted for "normal" margins of 2.54cm.

    for jj = shockvec                
        indSubplot = 1;
        for ii = varvec

            subplot(row,col,indSubplot);

            plot(steps,IRF{indRegimes}(:,ii,jj),'LineStyle','-','Color',rgb('dark blue'),'LineWidth',2);
            hold on

            plot(x_axis,'k','LineWidth',0.5)

            if ~isempty(bootIRFLowerCI)
                plot(steps,bootIRFLowerCI{indRegimes}(:,ii,jj),'LineStyle',':','Color',rgb('light blue'),'LineWidth',1.5);
            end
            hold on

            if ~isempty(bootIRFUpperCI)
                plot(steps,bootIRFUpperCI{indRegimes}(:,ii,jj),'LineStyle',':','Color',rgb('light blue'),'LineWidth',1.5);
            end

            xlim([1 numPer]);
            if ~isempty(vnames)
                title(vnames(ii,:), 'FontWeight','bold','fontsize',10); 
            end

            if any(ii == varvec(end-min(length(varvec),3)+1:end))
                xlabel('Quarters');
            end

            if any(ii == colvec)
                ylabel('Impulse response');
            end

            FigFont(fontsize);
            indSubplot = indSubplot + 1;
        end

        switch suptit
            case 1
                if ~isempty(vnames)
                    SupTitle( ['IRF - Shock to ' deblank(char(vnames(jj,:))) ' equation in Regime ' num2str(indRegimes)] );
                else
                    SupTitle( ['IRF - Shock to y_' num2str(jj) ' equation' ] );
                end
            case 0
                % Nothing
        end

        % Save
        set(gcf,'Color','w');
        FigName = [filename num2str(jj) num2str(indRegimes)];
        export_fig(FigName,'-pdf','-painters')
        clf('reset');
    end

    close all

end % for {indRegimes}


end % function {tvarxlirfplt}


% ---------------------------------
function [results,IRF] = checkInput(results,IRF)
% checkInput: Local function to check the validity of required inputs

% Define sizes and constants
numY = results(1).sizes.numY;

% results
validateattributes(results,{'struct'},{},'tvarxirfplt','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxirfplt:InputError','tvarxirfplt: The results structure provided must be returned by TVARX function')
end

% IRF
validateattributes(IRF,{'cell'},{},'tvarxlirfplt','IRF',1);	% check IRF
for indRegimes = 1:numel(IRF)
    validateattributes(IRF{indRegimes},{'numeric'},{'real','finite','size',[NaN numY numY]},'tvarxirfplt','IRF',1);
end


end % subfunction {checkInput}


% ---------------------------------
function [bootIRFLowerCI,bootIRFUpperCI,vnames,variables,shocks,filename,fontsize,suptit] = checkOptions(options,results,IRF) %#ok<STOUT>
% checkOptions: Local function to check the validity of options provided. Options are 
%               always optional. If the user does not provide any options, CHECKOPTIONS will
%               choose default ones.

% Define sizes and constants
numRegimes = numel(IRF);
numPer = cell(numRegimes,1);
numDims = cell(numRegimes,1);
numShocks = cell(numRegimes,1);
for indRegimes = 1:numRegimes
	[numPer{indRegimes},numDims{indRegimes},numShocks{indRegimes}] = size(IRF{indRegimes});
end
numExo  = results(1).sizes.numExo;
isCons  = results(1).constant;
isTrend = results(1).trend;

% Get default or user-provided options
getOptions(options, ...
    'bootIRFLowerCI',	cell.empty, ...       % IRF lower bound
	'bootIRFUpperCI',	cell.empty, ...       % IRF upper bound
    'vnames',           char.empty(numDims{1}+numExo-double(isCons)-double(isTrend),0), ...	% vector of variable names
    'variables',        true(1,numDims{1}), ...	% variables to be plotted
    'shocks',           true(1,numDims{1}), ...	% shocks to be plotted
    'filename',         'irf_', ...             % filename
    'fontsize',         16, ...                 % fontsize
    'suptit',           true);                  % superior title


% bootIRFLowerCI
switch ~isempty(bootIRFLowerCI)
    case 0 
        % Nothing
    case 1
        validateattributes(bootIRFLowerCI,{'cell'},{'numel',numRegimes},'tvarxirfplt','options.bootIRFLowerCI',3);
        for indRegimes = 1:numRegimes
            validateattributes(bootIRFLowerCI{indRegimes},{'numeric'},{'real','finite','size',[numPer{indRegimes} numDims{indRegimes} numShocks{indRegimes}]},'varxirfplt','options.bootIRFLowerCI',3);
        end
end

% bootIRFUpperCI
switch ~isempty(bootIRFUpperCI)
    case 0 
        % Nothing
    case 1
        validateattributes(bootIRFUpperCI,{'cell'},{'numel',numRegimes},'tvarxirfplt','options.bootIRFUpperCI',3);
        for indRegimes = 1:numRegimes
            validateattributes(bootIRFUpperCI{indRegimes},{'numeric'},{'real','finite','size',[numPer{indRegimes} numDims{indRegimes} numShocks{indRegimes}]},'varxirfplt','options.bootIRFUpperCI',3);
        end
end

% vnames
validateattributes(vnames,{'char'},{'nrows',numDims{1}+numExo-double(isCons)-double(isTrend)},'tvarxirfplt','options.vnames',3);

% variables
validateattributes(variables,{'logical'},{'size',[1 numDims{1}]},'tvarxirfplt','options.variables',3);

% shocks
validateattributes(shocks,{'logical'},{'size',[1 numDims{1}]},'tvarxirfplt','options.shocks',3);

% filename
validateattributes(filename,{'char'},{'nrows',1},'tvarxirfplt','options.filename',3);

% fontsize
validateattributes(fontsize,{'numeric'},{'real','finite','positive'},'tvarxirfplt','options.fontsize',3);

% suptit
validateattributes(suptit,{'logical'},{'numel',1},'tvarxirfplt','options.suptit',3);


end % subfunction {checkOptions}

