function results = mvrlsfit(y,x,Rest,r,algo)
% MVRLSFIT: Multivariate multiple linear regression using restricted (linear equality constrained) generalized least squares.
%
% Syntax: 
%
%       results = mvrlsfit(y,x,Rest,r,algo)
%
% Description:
%
%       Computes multivariate multiple linear regression using restricted (linear equality constrained) generalized least squares
%		of the form R*gamma - r = beta.
%
%       MVRLSFIT treats NaNs in X or Y as missing values, and removes them.
%
%       If columns of X are linearly dependent, select algo = {'qr','svd'}. MVRLSFIT sets the maximum possible number of elements of 
%       B to zero to obtain a "basic solution". In this case, the restriction [R,r] are not applied and the basic solution is returned.
%
% Input Arguments:
%
%       y                   -	[cell]      numRes-by-1 cell array, with each cell containing a (potentially) specific
%                                           numObs-by-numY matrix of dependent variable (or response observations, or 
%                                           regressands)
%       x                   -	[cell]      numRes-by-1 cell array, with each cell containing a numObs-by-numX 
%                                           independent variables matrix (regressors). Its a design matrix, 
%                                           with rows corresponding to observations and columns to predictor variables. 
%                                           Each 'x' cell must match the dimension (numObs) of the corresponding 'y' cell, 
%                                           the first 'x' being the design matrix for the first 'y', etc... 
%       Rest                -	[cell]      numRes-by-1 cell array, with each cell containing a restriction matrix 
%											((numX-by-numY)-by-((numX-by-numY)-numRest) where numRest < numX. numRest 
%											corresponds to the number of linear (equality) restrictions
%       r                   -	[cell]      numRes-by-1 cell array, with each cell containing a restriction vector ((numX-by-numY)-by-1)
%
% Optional Input Arguments:
%
%       algo                -   [char]      string, way to compute the OLS estimator of beta:
%                                               'explicit'  explicit formula for global minimum
%                                               'qr'        QR decomposition (default)
%                                               'svd'       SVD decomposition
%
% Output Arguments:
%
%       results, a structure with the following fields:
%           fitmeth     -   [char]      string, fitting method
%           description	-   [char]      string, method description
%           algo        -   [char]      string, solution algorithm
%           numRes      -	[integer]	scalar, number of response path
%           solmat      -   [struct]    structure, whose fields are matrices of the solution
%                                       If algo = 'explicit':
%                                       	nothing
%                                       If algo = 'qr':
%                                           Q           - 	[double]    (numObs-numNaN)-by-numX, Q from the QR decomposition of the design matrix
%                                           R           - 	[double]    numX-by-numX, R from the QR Decomposition of the design matrix    
%                                           perm        -   [double]    1-by-numX, permutation information
%                                       If algo = 'svd':
%                                           U           - 	[double]    (numObs-numNaN)-by-(numObs-numNaN), U from the SVD decomposition of the design matrix    
%                                           V           - 	[double]    numX-by-numX, V from the SVD decomposition of the design matrix
%                                           S           -   [double]    (numObs-numNaN)-by-numX, rectangular diagonal matrix with nonnegative real numbers on the diagonal
%           data        -	[struct]	structure, whose fields are matrices are the data
%                                       	y           -	[double]	numObs-by-numY, regressands
%                                           x           -	[double]	numObs-by-numX, regressors
%                                           ynan        -	[double]	(numObs-numNaN)-by-numY, regressands with NaN removed
%                                           xnan        -	[double]	(numObs-numNaN)-by-numX, regressors with NaN removed
%           beta        -	[double]	numX-by-numY, constrained (restricted) estimated coefficients
%           betaunc     -	[double]	numX-by-numY, unconstrained estimated coefficients
%           Rest        -	[double]	((numX-by-numY)-by-((numX-by-numY)-numRest), restriction matrix
%           r           -	[double]	((numX-by-numY)-by-1), restriction vector
%           sizes       -   [struct]    structure, whose fields are defined size
%                                           numObs      -	[integer]	scalar, number of observations
%                                           numObsNaN	-	[integer]	scalar, number of observations after NaN are removed
%                                           numY        -	[integer]	scalar, number of regressands
%                                           numX        -	[integer]	scalar, number of regressors
%                                           numNaN      -	[integer]	scalar, number of removed values (NaN)
%                                           numParams	-	[integer]	scalar, number of parameter
%                                           numActive	-	[integer]	scalar, number of active (unrestricted) parameters
%                                           numRest     -	[integer]	scalar, number of restrictions
%           nanvec      -	[logical]	vector of non-NaN logical indicator 
%           allNaN      -	[logical]	true if the series is all NaN
%           rankdef     -	[logical]	true if the series is rank deficient           
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       vec       
%
% References:
%
%       Roderick J. A. Little and Donald B. Rubin, Statistical Analysis with Missing Data, 2nd ed., John Wiley & Sons, Inc., 
%       2002.
%
%       Xiao-Li Meng and Donald B. Rubin, "Maximum Likelihood Estimation via the ECM Algorithm," Biometrika, Vol. 80, No. 2,
%       1993, pp. 267-278.
%
%       Joe Sexton and Anders Rygh Swensen, "ECM Algorithms that Converge at the Rate of EM," Biometrika, Vol. 87, No. 3,
%       2000, pp. 651-662.
%
%       A. P. Dempster, N.M. Laird, and D. B. Rubin, "Maximum Likelihood from Incomplete Data via the EM Algorithm," Journal 
%       of the Royal Statistical Society, Series B, Vol. 39, No. 1, 1977, pp. 1-37.
%
%       Lütkepohl, Helmut, "Econometric Analysis with Vector Autoregressive Models", EUI Working Papers, ECO 2007/11
%
% Notes:
%
%       none 
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2013
%           The generalized least squares section of the code is taken from IRIS toolbox by Jaromir Benes.


% Input and Output arguments checking
% ___________________________________

narginchk(4,5); nargoutchk(1,1);

validateattributes(y,   {'cell'},{},                'mvrlsfit','y',   1)
validateattributes(x,   {'cell'},{'numel',numel(y)},'mvrlsfit','x',   2)
validateattributes(Rest,{'cell'},{'numel',numel(y)},'mvrlsfit','Rest',3)
validateattributes(r,   {'cell'},{'numel',numel(y)},'mvrlsfit','r',   4)
switch nargin
    case 4
        algo = 'qr';
    case 5
        validateattributes(algo,{'char'},{},'rlsfit','algo',5)
end % switch {nargin} 


% Define size
% ___________

numRes = numel(y);


% Initialize the nested structures
% ________________________________

solmat = ...
    struct('Q',     [], ...         % (numObs-numNaN)-by-numX, Q from the QR decomposition of the design matrix
           'R',     [], ...         % numX-by-numX, R from the QR Decomposition of the design matrix    
           'perm',  [], ...         % 1-by-numX, permutation information
           'U',     [], ...         % (numObs-numNaN)-by-(numObs-numNaN), U from the SVD decomposition of the design matrix    
           'V',     [], ...         % numX-by-numX, V from the SVD decomposition of the design matrix
           'S',     []);            % (numObs-numNaN)-by-numX, rectangular diagonal matrix with nonnegative real numbers on the diagonal
data = ...
    struct('y',     [], ...         % numObs-by-1, regressand
           'x',     [], ...         % numObs-by-numX, regressors
           'ynan',	[], ...         % (numObs-numNaN)-by-1, regressand with NaN removed
           'xnan',	[]);            % (numObs-numNaN)-by-numX, regressors with NaN removed
sizes = ...
    struct('numObs',	[], ...     % scalar, number of observations before NaN are removed
           'numObsNaN',	[], ...     % scalar, number of observations after NaN are removed
           'numY',      [], ...     % scalar, number of regressands
           'numX',      [], ...     % scalar, number of regressors
           'numNaN',	[], ...     % scalar, number of removed values (NaN)
           'numParams',	[], ...     % scalar, number of parameters       
           'numActive',	[], ...     % scalar, number of active (unrestricted) parameters       
           'numRest',   []);        % scalar, number of restrictions

       
% Initialize the results structure with some fields that are constant across response paths
% _________________________________________________________________________________________

results(1:numRes,1) = ...
	struct('fitmeth',       'olsfit', ...	% string, fitting method
           'description',	'Multiple linear regression using ordinary least squares', ...	% string, method description
           'algo',          algo, ...   	% string, solution algorithm
           'numRes',        numRes, ...     % scalar, number of response path
           'solmat',        solmat, ...     % structure, whose fields are matrices of the solution
           'data',          data, ...       % structure, whose fields are matrices are the data
           'beta',          [], ...	        % numX-by-1, estimated coefficients
           'betaunc',       [], ...	        % numX-by-1, unconstrained estimated coefficients
           'Rest',          [], ...	        % numRest-by-numX, restriction matrix
           'r',             [], ...	        % numRest-by-1, restriction vector
           'sizes',         sizes, ...       % structure, whose fields are defined size
           'nanvec',        [], ...         % vector of non-NaN logical indicator 
           'allNaN',        false, ...      % logical, true if the series is all NaN
           'rankdef',       false);         % logical, true if the series is rank deficient           

       
% Loop over response vectors
% __________________________

for indRes = 1:numRes
    
    % Input and Output arguments checking
    validateattributes(y{indRes},   {'numeric'},{'real','nonempty','2d'},                                                     'mvrlsfit','y{indRes}',   1)
    validateattributes(x{indRes},   {'numeric'},{'real','nonempty','2d','nrows',size(y{indRes},1)},                           'mvrlsfit','x{indRes}',   2)
    validateattributes(Rest{indRes},{'numeric'},{'real','finite','nonempty','2d','nrows',size(x{indRes},2)*size(y{indRes},2)},'mvrlsfit','Rest{indRes}',3)
    validateattributes(r{indRes},   {'numeric'},{'real','finite','nonempty','2d','nrows',size(x{indRes},2)*size(y{indRes},2)},'mvrlsfit','r{indRes}',   4)

    % Define size
    [numObs,numY] = size(y{indRes});
    numX = size(x{indRes},2);
    numRest = size(Rest{indRes},1);

    % Remove missing values, if any
	ynan = y{indRes}(~any(isnan([y{indRes} x{indRes}]),2),:);
    xnan = x{indRes}(~any(isnan([y{indRes} x{indRes}]),2),:);
    nanvec = ~any(isnan([y{indRes} x{indRes}]),2);
    numNaN = numObs - length(ynan);
    numObsNaN = numObs - numNaN;

    % Restrictions
    RestRes = Rest{indRes};
    rRes = r{indRes};
    
    % Check if ynan is empty. If true, fill out the required fields and skip the regression for that column of 'y'
    if numObsNaN == 0
        results(indRes).data.y          = y{indRes};
        results(indRes).data.x          = x{indRes};
        results(indRes).data.ynan	    = ynan;
        results(indRes).data.xnan       = xnan;
        results(indRes).sizes.numObs	= numObs;
        results(indRes).sizes.numObsNaN = numObsNaN;
        results(indRes).sizes.numY	    = numY;
        results(indRes).sizes.numX	    = numX;
        results(indRes).sizes.numNaN    = numNaN;          
        results(indRes).nanvec          = nanvec;
        results(indRes).allNaN          = true;
        results(indRes).rankdef         = [];       
        continue
    end % if {numObsNaN == 0}
        
	% OLS estimation of beta
    switch lower(algo)
        case 'explicit' % Use explicit formula for global minimum
            betaunc = (xnan'*xnan)\xnan'*ynan;
        
        case 'qr'       % Use the rank-revealing QR decomposition to remove dependent columns of X, if any.
            [Q,R,perm] = qr(xnan,0);
            p = sum(abs(diag(R)) > max(numObsNaN,numX)*eps(R(1)));
            if p < numX
                warning('olsfit:RankDefDesignMat','olsfit: x is rank deficient to within machine precision');
                results(indRes).rankdef = true;
                R = R(1:p,1:p);
                Q = Q(:,1:p);
                perm = perm(1:p);
            end

            % Compute the OLS coefficients, filling in zeros in elements corresponding
            % to rows of X that were thrown out.
            betaunc = zeros(numX,numY);
            betaunc(perm,:) = R\(Q'*ynan);
            
        case 'svd'       % Use the SVD decomposition to identify dependent columns of X, if any.
            [U,S,V] = svd(xnan,0);
            s = diag(S);

            % Determine the effective rank p of A using singular values
            p = 1;
            while(p < size(xnan,2) && s(p+1) >= max(size(xnan))*eps*s(1))
                p = p+1;
            end
            if p < numX
                warning('olsfit:RankDefDesignMat','olsfit: x is rank deficient to within machine precision');
                results(indRes).rankdef = true;
            end
            
            d = 1./diag(S);
            betaunc = (V.*repmat(d',size(V,1),1))*(U'*ynan);

    end % switch {lower(algo)}
    
    % Compute the RLS coefficients
    if ~results(indRes).rankdef
        % Generalized least squares for parameter restrictions
        Omgi = eye(numY);
        beta = Inf;
        M = xnan.'*xnan;
        maxdiff = Inf;
        count = 0;
        while maxdiff > 1e-4 && count <= 100
            lastbeta = beta;
            c = vec(ynan') - kron(xnan,eye(numY))*rRes;
            % Estimate free hyperparameters.
            gamma = (RestRes.'*kron(M,Omgi)*RestRes)\(RestRes.'*kron(xnan.',Omgi)*c);
            % Compute parameters.
            beta = reshape(RestRes*gamma + rRes,[numY numX]);
            resid = ynan' - beta*xnan.';
            Omg = resid*resid.'/numObs;
            Omgi = inv(Omg);
            maxdiff = max(abs(beta(:) - lastbeta(:)));
            count = count + 1;
        end
        beta = beta.';
	else
        % If the rank is deficient, do not compute the restrictions
        beta = [];
        warning('mvrlsfit:RankDefDesignMat','mvrlsfit: x is rank deficient, constrained beta is not computed');
    end % if {~results(indRes).rankdef}

    % Put results inside results structure (help to keep the code readable)
	switch lower(algo)
        case 'explicit' % Use explicit formula for global minimum
            % Nothing
        
        case 'qr'       % Use the rank-revealing QR decomposition to remove dependent columns of X, if any.
            results(indRes).solmat.Q    = Q;
            results(indRes).solmat.R    = R;
            results(indRes).solmat.perm = perm;
        
        case 'svd'       % Use the SVD decomposition to identify dependent columns of X, if any.
            results(indRes).solmat.U = U;
            results(indRes).solmat.V = V;
            results(indRes).solmat.S = S;

	end % switch {lower(algo)}
	
    results(indRes).data.y          = y{indRes};
	results(indRes).data.x          = x{indRes};
	results(indRes).data.ynan       = ynan;
	results(indRes).data.xnan       = xnan;
	results(indRes).beta            = beta;
	results(indRes).betaunc         = betaunc;
	results(indRes).Rest            = RestRes;
	results(indRes).r               = rRes;
    results(indRes).sizes.numObs    = numObs;
    results(indRes).sizes.numObsNaN = numObsNaN;
    results(indRes).sizes.numY      = numY;
    results(indRes).sizes.numX      = numX;
	results(indRes).sizes.numNaN    = numNaN;  
	results(indRes).sizes.numParams	= numel(beta);  
	results(indRes).sizes.numActive	= numel(beta) - numRest;  
	results(indRes).sizes.numRest   = numRest;
	results(indRes).nanvec          = nanvec;
    
end % for {indRes}


end % function {mvrlsfit}

