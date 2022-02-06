function [bsdata,indices] = bootcont(data,lengthBlock,numBoot,numPer)
% BOOTCONT: Implements Politis' continuous-path block bootstrap for bootstrapping unit root series.
%
% Syntax: 
%
%       [bsdata,indices] = bootcont(data,numBoot,lengthBlock)
%
% Description:
%
%       Implements Politis'  continuous-path block bootstrap for bootstrapping unit root series.
%
%       To generate bootstrap sequences for other uses, such as bootstrapping vector processes, set data to (1:numObs)'.
%
% Input Arguments:
%
%       data        -	[double]	numObs-by-1 vector, data to be bootstrapped (should be unit root)
%       lengthBlock	-	[double]	scalar, average block length
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
%       E. Paparoditis, D.N.Politis "The Continuous-Path Block-Bootstrap", in Asymptotics in Statistics and Probability. 
%       Papers in honor of George Roussas, (Madan Puri, Ed.), VSP publications, Zeist (The Netherlands), 2001, pp. 305-320.
%
%       E. Paparoditis, D. N. Politis, "Residual-based Block Bootstrap for Unit Root Testing", Econometrica, Vol. 71, No. 3, pp. 813-855, 2003. 
%       (This is a revised version of the paper: `Unit Root Testing via the Continuous-Path Block Bootstrap', which is available as Discussion Paper 2001-06, 
%       Department of Economics, University of California, San Diego.) 
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

[numObs,numVars] = size(data);	% Get the length of the data
validateattributes(data,{'numeric'},{'real','finite','column'},'bootcont','data',1)
if numVars > 1
    error('bootcont:InputError','bootcont: data must be a column vector')
end
if numObs < 2 
    error('bootcont:InputError','bootcont: data must have at least 2 observations.')
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
numBlock = sign(data);
data = data.^2;

% Set up the bsdata and indices
bsdata = zeros(numPer,numBoot);
indices = zeros(numPer,numBoot);

% Initial positions
indices(1,:) = ceil(numObs*rand(1,numBoot));
bsdata(1,:) = data(indices(1,:))';

% Set up the random numbers
for i = 2:numPer
	select = rand(1,numBoot) < probBlock;
	indices(i,:) = ceil((numObs-1)*rand(1,numBoot)).*select+(mod(~select.*(indices(i-1,:)),numObs)+1);
	bsdata(i,:) = data(indices(i,:))';
end

for i = 1:numBoot
	scaling = 1;
	d = diff(indices);
	for j = 2:length(numBoot)-1;
        if d(j,i) ~= 1;
            scaling = bsdata(j,i)/bsdata(j-1,i);
        end
        bsdata(j,i) = bsdata(j,i)*scaling;
	end
end

% The indices make finding the bsdata simple
bsdata = sqrt(bsdata).*numBlock(indices);


end % function {bootcont}

