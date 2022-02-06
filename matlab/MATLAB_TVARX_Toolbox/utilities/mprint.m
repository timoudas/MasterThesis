function mprint(y,options)
% MPRINT: Print a matrix in formatted form.
%
% Syntax: 
%
%        mprint(y,options)
%
% Description:
%
%       Print time-series matrix or vector with dates and column labels.
%
%       Example:    options.cnames = char('col1','col2');
%                   options.rnames = char('rowlabel','row1','row2');
%                   mprint(y,options), prints entire matrix, column and row headings
%                   options2.endc = 3; options2.cnames = char('col1','col2','col3');
%                   or: mprint(y,options2), prints 3 columns of the matrix, just column headings 
%                   or: mprint(y), prints entire matrix, no column headings or row labels 
%
% Input Arguments:
%
%       y                   -   [double]	a numObs-by-numVars vector or matrix to be printed
%
% Optional Input Arguments:
%
%       options, a structure with the following possible fields and possible values:
%       options.endr       -	[integer]	scalar, ending row to print
%                                           (default = numObs)
%       options.begc       -	[integer]   scalar, beginning column to print
%                                           (default = 1)
%       options.endc       -	[integer]   scalar, ending column to print 
%                                           (default = numVars)        
%       options.cnames     -	[double]	a numVars-by-1 string vector of names for columns
%                                           ex: options.cnames = char('col1','col2');
%                                           (default = no column headings)
%       options.rnames     -   [char]       a (numObs+1)-by-1 string vector of names for rows
%                                           ex: options.rnames = char('Rows','row1','row2');
%                                           (default = no row labels)
%       options.fmt        -   [char]       a format string, 
%                                               ex: '%12.6f' or '%12d' 
%                                           (default = %10.4f)
%                                           or a numVars-by-1 string containing formats
%                                               ex: options.fmt = char('%12.6f','%12.2f','%12d'); for nvar = 3
%       options.fid        -   [char]       file-id for printing results to a file
%                                           (defaults to the MATLAB command window)
%                                               ex: fid = fopen('file.out','w'); 
%       options.rflag      -   [integer]	scalar, 1 for row #'s printed, 0 for no row #'s 
%                                           (default = 0) 
%       options.width      -   [integer]	scalar, number of columns before wrapping occurs 
%                                           (default = 80)        
%
% Output Arguments:
%
%       Printed matrix in formatted form
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
%       Defaults are used for options-elements not specified. Default wrapping occurs at 80 columns, which varies depending on the
%       format you use, e.g. %10.2f will wrap after 8 columns
%
% Copyright:
% 
%       (c) James P. LeSage, Dept of Economics
%           University of Toledo
%           2801 W. Bancroft St,
%           Toledo, OH 43606
%           jlesage@spatial-econometrics.com


% Input and Output arguments checking
% ___________________________________

narginchk(1,2); nargoutchk(0,0);

validateattributes(y,{'numeric'},{'2d','real'},'mprint','y',1);


% Options
% _______

fid = 1;
rflag = 0;
cflag = 0;
rnum = 0;
nfmts = 1;
cwidth = 80;
[nunObs nunVars] = size(y);
begr = 1;
endr = nunObs;
begc = 1;
endc = nunVars;
fmt = '%10.4f';

switch nargin
    case 1
        % rely on defaults
    case 2
        validateattributes(options,{'struct'},{},'mprint','options',2);
        fields = fieldnames(options);
        nf = length(fields);
        for i=1:nf
            if strcmp(fields{i},'fmt')
                fmts = options.fmt; 
                nfmts = size(fmts,1);
                if nfmts == nunVars
                    fmt = fmts;
                elseif nfmts == 1
                    fmt = fmts;
                else
                    error('mprint: wrong # of formats in string -- need nvar');
                end
            elseif strcmp(fields{i},'fid')
                fid = options.fid;
            elseif strcmp(fields{i},'begc');
                begc = options.begc;
            elseif strcmp(fields{i},'begr');
                begr = options.begr;
            elseif strcmp(fields{i},'endc');
                endc = options.endc;
            elseif strcmp(fields{i},'endr');
                endr = options.endr;
            elseif strcmp(fields{i},'width');
                cwidth = options.width;
            elseif strcmp(fields{i},'cnames');
                cnames = options.cnames;
                cflag = 1;
            elseif strcmp(fields{i},'rnames');
                rnames = options.rnames;
                rflag = 1;
            elseif strcmp(fields{i},'rflag');
                rnum = options.rflag;
            end
        end
end % switch {nargin}


% See if the user supplied row names and set rnum correct her mistake if she did this
if rflag == 1
    rnum = 0;
end


% Parse formats
% _____________

if nfmts == 1
	f1 = strtok(fmt,'%');
	f2 = strtok(f1,'.'); 
	if strcmp(f1,f2)
        f2 = strtok(f2,'d');
        dflag = 1;
        fflag = 0;
    else
        tmp1 = strtok(fmt,'f');
        tmp2 = strtok(tmp1,'.');
        tmp1 = tmp1(2:length(tmp1));
        tmp2 = tmp2(2:length(tmp2));
        opoint = num2str(str2num(tmp1) - str2num(tmp2)); %#ok<ST2NM>
        decimal = opoint(1,length(opoint));
        f2 = strtok(f2,'f');
        fflag = 1;
        dflag = 0;
	end
	f2 = str2num(f2); %#ok<ST2NM>
	nwide = floor(cwidth/f2); % 80 columns divided by format
	nvar = endc-begc+1;
	nsets = ceil(nvar/nwide);
else % wrapping in this case is based on widest format in the list
    nwidev = zeros(nfmts,1);
    nsetsv = zeros(nfmts,1);
    f2v = zeros(nfmts,1);
    dflagv = zeros(nfmts,1);
    fflagv = zeros(nfmts,1);
    decimalv = zeros(nfmts,1);
	for ii=1:nfmts;
        f1 = strtok(fmt(ii,:),'%');
        f2 = strtok(f1,'.');
        if strcmp(f1,f2)
            f2 = strtok(f2,'d');
            dflagv(ii,1) = 1;
            fflagv(ii,1) = 0;     
        else
            tmp1 = strtok(fmt(ii,:),'f');
            tmp2 = strtok(tmp1,'.');
            tmp1 = tmp1(2:length(tmp1));
            tmp2 = tmp2(2:length(tmp2));
            opoint = num2str(str2num(tmp1) - str2num(tmp2)); %#ok<ST2NM>
            decimalv(ii,1) = opoint(1,length(opoint));     
            f2 = strtok(f2,'f');
            fflagv(ii,1) = 1;
            dflagv(ii,1) = 0;     
        end
        f2v(ii,1) = str2num(f2); %#ok<ST2NM>
        nwidev(ii,1) = floor(cwidth/f2v(ii,1)); % cwidth columns divided by format
        nvar = endc-begc+1;
        nsetsv(ii,1) = ceil(nvar/nwidev(ii,1));   
	end
    nsets = min(nsetsv); 
    nwide = max(nwidev);
end

% If we have row and column labels adjust variable labels and column heading strings
% to match the width of the printing format
if rnum == 1
    dstr = 'Obs#';
end

if cflag == 1 % we have column headings
	[vsize nsize] = size(cnames); % error check cnames argument
	if vsize ~= nunVars
        error('mprint:InputError','mprint: Wrong number of cnames in mprint');
	end   
	if nfmts == 1 % case of only 1 format string
        nmax = max(f2,nsize); % build format strings, % based on widest format 
        sfmt = ['%', num2str(nmax)];
        sfmt = [sfmt,'s ']; 
        ffmt = ['%', num2str(nmax)];
        if dflag == 1
            ffmt = [ffmt,'d '];
        elseif fflag == 1
            ffmt = [ffmt,'.'];
            ffmt = [ffmt,decimal];
            ffmt = [ffmt,'f '];
        end
    else % we have multiple format strings, process each
        sfmtv = []; 
        for ii = 1:nfmts % find and parse multiple formats
            nmax = max(f2v(ii,:),nsize); % build format strings, based on widest format 
            sfmtv{ii} = ['%', num2str(nmax)]; %#ok<AGROW>
            sfmtv{ii} = [sfmtv{ii},'s '];  %#ok<AGROW>
            ffmtv{ii} = ['%', num2str(nmax)]; %#ok<AGROW>
            if dflagv(ii,1) == 1
                ffmtv{ii} = [ffmtv{ii},'d ']; %#ok<AGROW>
            elseif fflagv(ii,1) == 1
                ffmtv{ii} = [ffmtv{ii},'.']; %#ok<AGROW>
                ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];  %#ok<AGROW>   
                ffmtv{ii} = [ffmtv{ii},'f ']; %#ok<AGROW>
            end
        end % end of for ii loop
	end % end of if-else
elseif cflag == 0 % we have no column headings
	if nfmts == 1 % case of only 1 format string
        nmax = f2; % augment format string with a space (the hard way) 
        ffmt = ['%', num2str(nmax)];
        if dflag == 1
            ffmt = [ffmt,'d '];
        elseif fflag == 1
            ffmt = [ffmt,'.'];
            ffmt = [ffmt,decimal];
            ffmt = [ffmt,'f '];
        end
    else % we have multiple format strings, process each
        sfmtv = []; 
        for ii = 1:nfmts % find and parse multiple formats
            nmax = f2v(ii,:); % augment format strings with a space 
            ffmtv{ii} = ['%', num2str(nmax)]; %#ok<AGROW>
            if dflagv(ii,1) == 1
                ffmtv{ii} = [ffmtv{ii},'d ']; %#ok<AGROW>
            elseif fflagv(ii,1) == 1
                ffmtv{ii} = [ffmtv{ii},'.']; %#ok<AGROW>
                ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];   %#ok<AGROW>  
                ffmtv{ii} = [ffmtv{ii},'f ']; %#ok<AGROW>
            end
        end % end of for ii loop
	end % end of if-else    
end % end of if-elseif cflag == 0,1
   
if rflag == 1 % we have row labels
	[vsize nsize] = size(rnames); % error check cnames argument
	if vsize ~= nunObs+1
        error('mprint:InputError','mprint: Wrong number of rnames in mprint');
	end 
    rfmt = ['%', num2str(nsize)]; 
    rfmt = [rfmt,'s ']; 
end % end of if rflag == 1

if rflag == 0 && cflag == 0
    ffmt = fmt;
end


% Print matrix
% ____________

for j = 1:nsets
	if nfmts == 1 % print row header and column headers
        if rnum == 1
            fprintf(fid,'%5s ',dstr);     
        elseif rflag == 1    
            fprintf(fid,rfmt,rnames(1,:));
        end  
        if cflag == 1
            for i = (j-1)*nwide+begc:j*nwide+begc-1
                if i <= endc
                    % find version #; 
                    %[version,junk] = version; vers = str2num(version);
                    %if vers == 5.2
                        fprintf(fid,sfmt,strjust(cnames(i,:),'right'));
                    %else
                        %fprintf(fid,sfmt,strjust(cnames(i,:)));
                    %end
                end
            end
        end
        fprintf(fid,'\n');
	else % we have multiple formats
        if rnum == 1
            fprintf(fid,'%5s ',dstr);     
        elseif rflag == 1   
            fprintf(fid,rfmt,rnames(1,:));
        end
        if cflag == 1
            for i = (j-1)*nwide+begc:j*nwide+begc-1
                if i <= endc
                    % find version #; 
                    %[version,junk] = version; vers = str2num(version);
                    %if vers == 5.2
                        fprintf(fid,sfmtv{i},strjust(cnames(i,:),'right'));
                    %else
                        %fprintf(fid,sfmtv{i},strjust(cnames(i,:)));
                    %end
                end
            end
        end
        fprintf(fid,'\n');
	end % end of if-else nfmts
	for k = begr:endr % print row labels and numbers in matrix
        if rnum == 1
            fprintf(fid,'%5d ',k);
        elseif rflag == 1        
            fprintf(fid,rfmt,rnames(k+1,:));
        end
        for l = (j-1)*nwide+begc:j*nwide+begc-1
            if l <= endc
                if nfmts == 1
                    fprintf(fid,ffmt,y(k,l));
                else
                    fprintf(fid,ffmtv{l},y(k,l));
                end
            end
        end % end of for l
        fprintf(fid,'\n');
	end % end of for k
%    fprintf(fid,'\n');
end % end of for j


end % function {mprint}

