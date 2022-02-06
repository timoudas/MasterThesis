function bsdata = bootiid(data,numBoot,numPer)
% BOOTIID: Implements an iid bootstrap resample.
%
% Syntax: 
%
%       bsdata = bootiid(data,numBoot)
%
% Description:
%
%       Implements bootstrap, iid, resample. To generate bootstrap sequences for other uses, such as bootstrapping 
%       vector processes, simple set data to (1:1:numObs)'. For a vector input data of size [numObs,1], the  resampling 
%       procedure produces a matrix of size [numObs,numBoot] with columns being resamples of the input vector.
%
%       For a matrix input data of size [numObs,numVars], the resampling procedure produces a 3D matrix of size 
%       [numObs,numVars,numBoot] with bsdata(:,:,i), i = 1,...,numBoot, being a resample of the input matrix.
%
% Input Arguments:
%
%       data        -	[double]	numObs-by-1 vector, data to be bootstrapped
%
%       where
%           numObs is the sample size
%
% Optional Input Arguments:
%
%       numBoot     -	[integer]	scalar, number of bootstraps (default = 999)
%       numPer      -	[integer]	scalar, number of period to generate (default = numObs)
%
% Output Arguments:
%
%       bsdata      -	[double]	numPer-by-numBoot array, bootstrapped data
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
%       Efron, B.and Tibshirani, R.  An Introduction to the Bootstrap. Chapman and Hall, 1993.
%
%       Zoubir, A.M. Bootstrap: Theory and Applications. Proceedings of the SPIE 1993 Conference on Advanced 
%       Signal Processing Algorithms, Architectures and Implementations. pp. 216-235, San Diego, July  1993.
%
%       Zoubir, A.M. and Boashash, B. The Bootstrap and Its Application in Signal Processing. IEEE Signal 
%       Processing Magazine, Vol. 15, No. 1, pp. 55-76, 1998.
%
% Notes:
%
%       none
%
% Copyright:
% 
%       (c) Created by A. M. Zoubir and D. R. Iskander, May 1998


% Input and Output arguments checking
% ___________________________________

narginchk(1,3); nargoutchk(0,1);

[numObs,numVars] = size(data);	% Get the length of the data

validateattributes(data,{'numeric'},{'real','finite','column'},'bootiid','data',1)
if numVars > 1
    error('bootiid:InputError','bootiid: data must be a column vector')
end
if numObs < 2 
    error('bootiid:InputError','bootiid: data must have at least 2 observations.')
end


% Options
% _______

switch nargin
    case 1
        numBoot = 999;
        numPer = numObs;
    case 2
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootiid','numBoot',2)
        numPer = numObs;
    case 3
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootiid','numBoot',2)
        validateattributes(numPer, {'numeric'},{'real','finite','positive','scalar','integer'},'bootiid','numPer', 3)
end


% Bootstrap
% _________

bsdata = data(ceil(numObs*rand(numPer,numBoot)));    


end % function {bootiid}



