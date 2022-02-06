function [isStable,eigARmax] = tvarxstab(results)
% TVARXSTAB: Test TVARX model for stability.
%
% Syntax:
%
%       [isStable,eigAR] = tvarxstab(results)
%
% Description:
%
%       Test TVARX model for stability.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by TVARX function
%
% Optional Input Arguments:
%
%       none       
%
% Output Arguments:
%
%       isStable    -	[logical]   Boolean vector indicating the stability of tvarx model. A value of true indicates stability and 
%                                   a value of false indicates instability. One indicator by regime.
%       eigAR       -	[double]    vector, eigenvalue of largest magnitude for the AR component of results. The process is AR-stable if 
%                                   and only if eigAR < 1. One eigen value by regime.
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, tvarxcompmat
%
% References:
%
%       none
%
% Notes:
%
%       A process with nonconstant exogenous inputs can be unstable or singular. This cannot be determined from results unless there are no
%       exogenous inputs. TVARXSTAB determines if the AR components of a TVARX model are stable without considering exogenous inputs. 
%       In this sense, it may be more appropriate to speak of "AR-stability" (if the AR component is stable).
%
%       A stable TVARX process is stationary but the converse is not true. Although the terms "stable," "stationary," and "covariance-
%       stationary" are often used interchangeably, a process is stationary if and only if its first and second moments are independent of 
%       time.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2013


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,2);

results = checkInput(results);


% Test stability
% ______________

CompMat = tvarxcompmat(results);
eigARmax = zeros(size(results,1),1);
for indRegimes = 1:size(results,1)
    eigARvalues = abs(eig(CompMat{indRegimes}));
    eigARmax(indRegimes) = max(eigARvalues);
end


% Determine if the processes are stable
% _____________________________________

isStable = (eigARmax < 1);


end % function {tvarxstab}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'tvarxstab','results',1);
if ~strcmpi(results(1).meth,'tvarx')
    error('tvarxstab:InputError','tvarxstab: The results structure provided must be returned by TVARX function')
end

end % subfunction {checkInput}

