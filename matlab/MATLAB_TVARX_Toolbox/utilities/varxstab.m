function [isStable,eigARmax] = varxstab(results)
% VARXSTAB: Test VARX model for stability.
%
% Syntax:
%
%       [isStable,eigAR] = varxstab(results)
%
% Description:
%
%       Test VARX model for stability.
%                
% Input Arguments:
%
%       results     -   [struct]    a results structure returned by VARX function
%
% Optional Input Arguments:
%
%       none       
%
% Output Arguments:
%
%       isStable    -	[logical]   Boolean scalar indicating the stability of varx model. A value of true indicates stability and 
%                                   a value of false indicates instability.
%       eigAR       -	[double]    scalar, eigenvalue of largest magnitude for the AR component of results. The process is AR-stable if 
%                                   and only if eigAR < 1.
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       checkInput, varxcompmat
%
% References:
%
%       none
%
% Notes:
%
%       A process with nonconstant exogenous inputs can be unstable or singular. This cannot be determined from results unless there are no
%       exogenous inputs. VARXSTAB determines if the AR components of a VARX model are stable without considering exogenous inputs. 
%       In this sense, it may be more appropriate to speak of "AR-stability" (if the AR component is stable).
%
%       A stable VARX process is stationary but the converse is not true. Although the terms "stable," "stationary," and "covariance-
%       stationary" are often used interchangeably, a process is stationary if and only if its first and second moments are independent of 
%       time.
%
% Copyright:
% 
%       (c) Gabriel Bruneau, 2012-2014


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,2);

results = checkInput(results);


% Test stability
% ______________

CompMat = varxcompmat(results);
eigARvalues = abs(eig(CompMat));
eigARmax = max(eigARvalues);

% Determine if the processes are stable
isStable = (eigARmax < 1); 


end % function {varxstab}


% ---------------------------------
function results = checkInput(results)
% checkInput: Local function to check the validity of required inputs

% results
validateattributes(results,{'struct'},{},'varxstab','results',1);
if ~strcmpi(results.meth,'varx')
    error('varxstab:InputError','varxstab: The results structure provided must be returned by VARX function')
end

end % subfunction {checkInput}


