function maxfig(hFig)
% MAXFIG: Maximize a figure window to fill the entire screen.
%
% Syntax: 
%
%       maximize
%       maximize(hFig)
%
% Description:
%
%       Maximizes the current or input figure so that it fills the whole of the
%       screen that the figure is currently on. This function is platform
%       independent.
%
% Input Arguments:
%
%       none
%
% Optional Input Arguments:
%
%       hFig	-	[figure]    Handle of figure to maximize. 
%                                   Default: gcf (get current figure)
%
% Output Arguments:
%
%       Maximizes the current or input figure.
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
%       (c) Oliver Woodford, 2009


% Input and Output arguments checking
% ___________________________________

narginchk(0,1); nargoutchk(0,0);

% Options
% _______

if nargin < 1
    hFig = gcf;
end

% Resize figure
% _____________

drawnow % Required to avoid Java errors

jFig = get(handle(hFig), 'JavaFrame'); 
jFig.setMaximized(true);


end % function {maxfig}

