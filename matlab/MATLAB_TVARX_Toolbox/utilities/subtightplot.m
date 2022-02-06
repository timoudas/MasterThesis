function h = subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% SUBTIGHTPLOT: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
%               neighbouring subplots.
%
% Syntax: 
%
%        h = subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%
% Description:
%
%       A wrapper function for Matlab function subplot. Adds the ability to define the gap between
%       neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
%       subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input Arguments:
%
%       m,n,p	-   [integer]	breaks the Figure window into an m-by-n matrix of small axes, selects the p-th axes for
%                               the current plot, and returns the axes handle.  The axes are
%                               counted along the top row of the Figure window, then the second
%                               row, etc.
%
% Optional Input Arguments:
%
%       gap     -	[double]    two elements vector [vertical,horizontal] defining the gap between neighbouring axes. 
%                               Default : 0.01. (Note this vale will cause titles legends and labels to collide with the subplots, 
%                               while presenting relatively large axis.) 
%       marg_h  -	[double]    margins in height in normalized units (0...1) or [lower uppper] for different lower and upper margins 
%       marg_w  -	[double]    margins in width in normalized units (0...1) or [left right] for different left and right margins 
%
% Output Arguments:
%
%       Same as subplot:
%       none, or axes handle according to function call
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       getOptions
%
% References:
%
%       none
%
% Notes:
%
%       Note that if additional elements are used in order to be passed to subplot, gap parameter must
%       be defined. For default gap value use empty element - [].
%
% Copyright:
% 
%       (c)	Felipe G. Nievinski, 2012
%           Pekka Kumpulainen, 2010
%           Nikolay S., 2011


% Input and Output arguments checking
% ___________________________________

narginchk(3,inf); nargoutchk(0,1);

% Note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row] = ind2sub([n,m],p);  

% Note subplot suppors vector p inputs - so a merged subplot of higher dimentions will be created
subplot_cols = 1 + max(subplot_col) - min(subplot_col); % number of column elements in merged subplot 
subplot_rows = 1 + max(subplot_row) - min(subplot_row); % number of row elements in merged subplot   


% Options
% _______

if (nargin < 4) || isempty(gap)
    gap = 0.01;
end
if (nargin < 5) || isempty(marg_h)
    marg_h = 0.05;
end
if (nargin < 5) || isempty(marg_w)
    marg_w = marg_h;
end
if isscalar(gap)
    gap(2) = gap;
end
if isscalar(marg_h)
    marg_h(2) = marg_h;
end
if isscalar(marg_w)
    marg_w(2) = marg_w;
end
gap_vert   = gap(1);
gap_horz   = gap(2);
marg_lower = marg_h(1);
marg_upper = marg_h(2);
marg_left  = marg_w(1);
marg_right = marg_w(2);


% Subplot
% _______

% Single subplot dimensions
%height = (1-(m+1)*gap_vert)/m;
%axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
height = (1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
%width = (1-(n+1)*gap_horz)/n;
%axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
width = (1-(marg_left+marg_right)-(n-1)*gap_horz)/n;

% Merged subplot dimensions:
merged_height = subplot_rows*(height+gap_vert) - gap_vert;
merged_width = subplot_cols*(width +gap_horz) - gap_horz;

% Merged subplot position:
merged_bottom = (m-max(subplot_row))*(height+gap_vert) +marg_lower;
merged_left = (min(subplot_col)-1)*(width+gap_horz) +marg_left;
pos_vec = [merged_left merged_bottom merged_width merged_height];

% h_subplot = subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h = subplot('Position',pos_vec,varargin{:});

if nargout < 1
    clear h
end


end % function  {subtightplot}

