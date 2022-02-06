function [bsdata,indices] = bootstationary(data,lengthBlock,numBoot,numPer)
% BOOTSTATIONARY: Implements a stationary bootstrap for bootstrapping stationary, dependent series.
%
% Syntax: 
%
%       [bsdata,indices] = bootstationary(data,numBoot,lengthBlock)
%
% Description:
%
%       Implements a stationary bootstrap for bootstrapping stationary, dependent series.
%
%       To generate bootstrap sequences for other uses, such as bootstrapping vector processes, set data to (1:numObs)'.
%
% Input Arguments:
%
%       data        -	[double]	numObs by 1 vector, data to be bootstrapped
%       lengthBlock	-	[double]	scalar, average block length. probBlock, the probability of starting a new block, is 
%                                   defined probBlock = 1/lengthBlock
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
%       bsdata      -	[double]	numPer-by-numBoot matrix, bootstrapped data
%       indices     -	[double]	numPer-by-numBoot matrix, locations of the original bsdata = data(indices);
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
%       D.N.Politis, J.P.Romano, "The Stationary Bootstrap", J. Amer. Statist. Assoc., vol. 89, No. 428, December 1994, pp. 1303-1313.
%
%       D.N.Politis, J.P.Romano, "Limit Theorems for Weakly Dependent Hilbert Space Valued Random Variables with Applications to the Stationary Bootstrap", 
%       Statistica Sinica, vol. 4, No. 2, July 1994, pp. 461-476. 
%
% Notes:
%
%       none
%
% Copyright:
% 
%       (c) Kevin Sheppard, 2001 (kevin.sheppard@economics.ox.ac.uk)


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(0,2);

[numObs numVars] = size(data);	% Get the length of the data

validateattributes(data,{'numeric'},{'real','finite','column'},'bootstationary','data',1)
if numVars > 1
    error('bootstationary:InputError','bootstationary: data must be a column vector')
end
if numObs < 2 
    error('bootstationary:InputError','bootstationary: data must have at least 2 observations.')
end
validateattributes(lengthBlock,{'numeric'},{'real','finite','positive','scalar','>=',1,'<=',numObs},'bootstationary','lengthBlock',2)


% Options
% _______

switch nargin
    case 1
        numBoot = 999;
        numPer = numObs;
    case 2
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootstationary','numBoot',3)
        numPer = numObs;
    case 3
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootstationary','numBoot',3)
        validateattributes(numPer, {'numeric'},{'real','finite','positive','scalar','integer'},'bootstationary','numPer', 4)
end


% Bootstrap
% _________

% Define the probability of a new block
probBlock = 1/lengthBlock;

% Set up the bsdata and indices
indices = zeros(numPer,numBoot);

% Initial positions
indices(1,:) = ceil(numObs*rand(1,numBoot));

% Set up the random numbers
select = rand(numPer,numBoot) < probBlock;
indices(select) = ceil(rand(1,sum(sum(select)))*numObs);
for i = 2:numPer
    % Determine whether we stay (rand > probBlock) or move to a new starting value (rand < probBlock)
    indices(i,~select(i,:)) = indices(i-1,~select(i,:)) + 1;
end

% The indices make finding the bsdata simple
if any(any(indices > numObs))
    check = 0;
else
    check = 1;
end

while check == 0
	indices(indices > numObs) = indices(indices > numObs) - numObs;
    if any(any(indices > numObs))
        check = 0;
    else
        check = 1;
    end
end
bsdata = data(indices);


end % function {bootstationary}

