function [bsdata,indices] = bootblock(data,lengthBlock,numBoot,numPer)
% BOOTBLOCK: Implements a circular block bootstrap for bootstrapping stationary, dependent series.
%
% Syntax: 
%
%       [bsdata,indices] = bootblock(data,lengthBlock,numBoot,numPer)
%
% Description:
%
%       Implements a circular block bootstrap for bootstrapping stationary, dependent series.
%
%       To generate bootstrap sequences for other uses, such as bootstrapping vector processes, set data to (1:numObs)'.
%
% Input Arguments:
%
%       data        -	[double]	numObs-by-1 vector, data to be bootstrapped
%       lengthBlock	-	[integer]	scalar, block length
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
%       bsdata      -	[double]	numObs-by-numBoot matrix, bootstrapped data
%       indices     -	[double]	numObs-by-numBoot matrix, locations of the original bsdata = data(indices);
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
%       (c) Kevin Sheppard, 2001 (kevin.sheppard@economics.ox.ac.uk)


% Input and Output arguments checking
% ___________________________________

narginchk(2,4); nargoutchk(0,2);

[numObs numVars] = size(data);	% Get the length of the data
validateattributes(data,{'numeric'},{'real','finite','column'},'bootblock','data',1)
if numVars > 1
    error('bootblock:InputError','bootblock: data must be a column vector')
end
if numObs < 2 
    error('bootblock:InputError','bootblock: data must have at least 2 observations.')
end
validateattributes(lengthBlock,{'numeric'},{'real','finite','positive','scalar','integer','>=',1,'<=',numObs},'bootblock','lengthBlock',2)


% Options
% _______

switch nargin
    case 1
        numBoot = 999;
        numPer = numObs;
    case 2
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootblock','numBoot',3)
        numPer = numObs;
    case 3
        validateattributes(numBoot,{'numeric'},{'real','finite','positive','scalar','integer'},'bootblock','numBoot',3)
        validateattributes(numPer, {'numeric'},{'real','finite','positive','scalar','integer'},'bootblock','numPer', 4)
end


% Bootstrap
% _________

% Make it easy to do the circular bootstrap
data = [data;data(1:lengthBlock-1)];

% Compute the number of blocks needed
%numBlock = ceil(numObs/lengthBlock);
numBlock = ceil(numPer/lengthBlock);

% Generate the starting points
Bs = floor(rand(numBlock,numBoot)*numObs)+1;
indices = zeros(numBlock*lengthBlock,numBoot);
index = 1;

% Adder is a variable that needs to be added each loop
adder = repmat((0:lengthBlock-1)',1,numBoot);
%for i = 1:lengthBlock:numObs
for i = 1:lengthBlock:numPer
    indices(i:(i+lengthBlock-1),:) = repmat(Bs(index,:),lengthBlock,1) + adder;
    index = index+1;
end

% The indices make finding the bsdata simple
indices = indices(1:numPer,:);
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


end % function {bootblock}
    
