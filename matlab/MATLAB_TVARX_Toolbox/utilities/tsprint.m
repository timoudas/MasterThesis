function tsprint(y,cstruct,options)
% TSPRINT: Print time-series matrix or vector with dates and column labels.
%
% Syntax: 
%
%        tsprint(y,cstruct), entire series, no variable names, or
%        tsprint(y,cstruct,options)
%
% Description:
%
%       Print time-series matrix or vector with dates and column labels.
%
%       Example:    cstr = cal(1980,1,12);
%                   tsprint(y,cstr,struct('begp',13,'endp',24), would print data for 1981
%                   tsprint(y,cstr,struct('begp',ical(1981,1,cstr),'endp',ical(1981,12,cstr)), which would print the same data for 1981
%
% Input Arguments:
%
%       y       -   [double]	matrix (or vector) of series to be printed
%       cstruct	-   [struct]	structure, returned by cal()
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       begp	-   [integer]	scalar, the beginning observation to print
%       endp 	-   [integer]	scalar, the ending period to print
%       vnames	-   [char]      a string matrix of names for a header
%                                   ex: vnames = char('y','x1','x2','x3');
%       fmt     -   [char]      a format string
%                                   ex: '%12.6f' or '%12d'
%
% Output Arguments:
%
%       Printed time-series matrix or vector with dates and column labels on the command window
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       getOptions, mprint
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
%       (c) James P. LeSage, Dept of Economics
%           University of Toledo
%           2801 W. Bancroft St,
%           Toledo, OH 43606
%           jlesage@spatial-econometrics.com


% Input and Output arguments checking
narginchk(2,3); nargoutchk(0,0);

validateattributes(y,{'numeric'},{'2d','real'},'tsprint','y',1);
validateattributes(cstruct,{'struct'},{},'tsprint','cstruc',2);

% Options
if nargin < 3
    options = struct();
elseif nargin == 3
	validateattributes(options,{'struct'},{},'tsprint','options',3);
end
getOptions(options, ...
    'begp',     1, ...
    'endp',     size(y,1), ...
    'vnames',   [], ...
    'fmt',      '%10.4f');     

validateattributes(begp,{'numeric'},{'real','finite','scalar','>=',1},'tsprint','options.begp',1); 
validateattributes(endp,{'numeric'},{'real','finite','scalar','<=',size(y,1)},'tsprint','options.endp',1); 
nflag = 0;
if ~isempty(vnames) %#ok<NODEF>
    nflag = 1;
    validateattributes(vnames,{'char'},{'nrows',size(y,2)},'tsprint','options.vnames',3);
end
validateattributes(fmt,{'char'},{},'tsprint','options.fmt',3); 

% Variables names
numVars = size(y,2);
if nflag == 0       % no variable names supplied, make some up
	vnames = [];
	for i = 1:numVars
        snames = 'series';
        name = [snames num2str(i)];
        vnames = char(vnames,name);
	end
    vnames = vnames(2:end,:);   % remove the empty space on the first row of vnames
else
	vsize = size(vnames,1);     % error check vnames argument
	if vsize ~= numVars
        error('tsprint:InputError','tsprint: wrong number of vnames in tsprint'); 
	end
end

% Dates
rnames = 'Date';
for k = begp:endp
    rnames = char(rnames,tsdate(cstruct,k));
end

% Print timeseries
in.rnames = rnames;
in.fmt = fmt;
in.cnames = vnames;
mprint(y(begp:endp,:),in);

end % function {tsprint}
 
