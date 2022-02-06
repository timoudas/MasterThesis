function results = mvlinreg(y,x,method,algo,varargin)
% MVLINREG: Multivariate multiple linear regression.
%
% Syntax: 
%
%       results = mvlinreg(y,x,'ols',algo)
%       results = mvlinreg(y,x,'rls',Rest,r,algo)
%
% Description:
%
%       Computes multivariate multiple linear regression. MVLINREG treats NaNs in X or Y as missing values, and removes them.
%
% Input Arguments:
%
%       y                   -	[cell]      numRes-by-1 cell array, with each cell containing a (potentially) specific
%                                           numObs-by-numY matrix of dependent variable (or response observations, or 
%                                           regressands)
%       x                   -	[cell]      numRes-by-1 cell array, with each cell containing a numObs-by-numVars 
%                                           independent variables matrix (regressors). Its a design matrix, 
%                                           with rows corresponding to observations and columns to predictor variables. 
%                                           Each 'x' cell must match the dimension (numObs) of the corresponding 'y' cell, 
%                                           the first 'x' being the design matrix for the first 'y', etc... 
%       method              -   [char]      string, estimation method
%                                               'ols'	ordinary least squares
%                                               'rls'   restricted (equality constrained) generalized least squares
%
% Optional Input Arguments:
%
%       algo                -   [char]      string, way to compute the OLS estimator of beta:
%                                               'explicit'  explicit formula for global minimum
%                                               'qr'        QR decomposition
%                                               'svd'       SVD decomposition
%       varargin, variable arguments in for the specific called function
%       When method = 'ols', no varargin
%       When method = 'rls', in this order:
%           Rest            -	[cell]      numRes-by-1 cell array, with each cell containing a restriction matrix (numRest x numX) 
%                                           where numRest < numX. numRest corresponds to the number of linear (equality) restrictions
%           r               -	[cell]      numRes-by-1 cell array, with each cell containing a restriction vector (numRest x 1)
%
% Output Arguments:
%
%       results, a structure specific to the estimation method used          
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
%       This is simply a wrapper function that calls another function. 
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2013


% Input and Output arguments checking
% ___________________________________

narginchk(3,Inf); nargoutchk(1,1);

validateattributes(method,  {'char'},{},'mvlinreg','method',3)
if nargin >= 4
    validateattributes(algo,{'char'},{},'mvlinreg','algo',  4)
else
    algo = 'qr';
end


% Estimation
% __________

switch method
    case 'ols'
        results = mvolsfit(y,x,algo);
    case 'rls'
        results = mvrlsfit(y,x,varargin{1},varargin{2},algo);
    otherwise
        error('mvlinreg:InputError','mvlinreg: method algorithm must be ''ols'' or ''rls''')
end % switch {method}


end % function {mvlinreg}

