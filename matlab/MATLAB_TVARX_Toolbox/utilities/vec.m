function datavec = vec(data,appendDirection)
% VECC: Creates a column vector by appending the rows or columns of a matrix to each other.
%
% Syntax: 
%
%       datavec = vec(data), or
%       datavec = vec(data,'rows')
%
% Description:
%
%       Creates a column vector by appending the rows or columns of a matrix to each other.
%
%       Example:    x = [1 2;
%                        3 4];
%                   yc = vec(x);
% 
%                   x = 1.000000 2.000000 
%                       3.000000 4.000000 
%                   yc = 1.000000 
%                        3.000000 
%                        2.000000 
%                        4.000000
%                   
%                   x = [1 2;
%                        3 4];
%                   yc = vec(x,'rows');
% 
%                   x = 1.000000 2.000000 
%                       3.000000 4.000000 
%                   yc = 1.000000 
%                        2.000000 
%                        3.000000 
%                        4.000000
%
% Input Arguments:
%
%       data	-	[double]    row-by-col, input matrix or vector
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       datavec	-	[double]    row*col-by-1 vector or matrix, the rows or columns of data appended to each other
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
%       TO DO: 1) Fix the 'rows' case when data is a multidimensional array (more than 2-D)
%
% Copyright:
% 
%       (c) Gabriel Bruneau (2011)


% Input and Output arguments checking
narginchk(1,2); nargoutchk(0,1);

% Apply the operator vec to rows or columns
switch nargin
    case 1
        vecDirection = 'cols';
    case 2
        vecDirection = appendDirection;
end % switch {nargin}       

% vec operator
switch vecDirection
    case 'cols' % append by columns
        datavec = data(:);

        % or (same output, but slower)
        %datavec = reshape(data,size(data,1)*size(data,2),1);

    case 'rows' % append by rows
        
        % Check if data is vector or matrix
        if ~ismatrix(data)
            error('vec:InputError','vec: data must be a vector or matrix')
        end

        datatemp = data.';
        datavec = vec(datatemp); 

        % or (same output, but slower)
        %datavec = reshape(data.',size(data,1)*size(data,2),1);        
end % switch {vecDirection}
        
end % function {vec}

