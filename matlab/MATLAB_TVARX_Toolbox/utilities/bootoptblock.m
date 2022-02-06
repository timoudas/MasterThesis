function Bstar = bootoptblock(data)
% BOOTOPTBLOCK: Optimal block length for the stationary bootstrap or circular block bootstrap.
%
% Syntax: 
%
%       Bstar = bootoptblock(data)
%
% Description:
%
%       This is a function to select the optimal (in the sense of minimising the MSE of the estimator of the long-run
%       variance) block length for the stationary bootstrap or circular block bootstrap.
%
%       The optimal average block length for the stationary bootstrap, and it does not need to be an integer. 
%       The optimal block length for the circular block bootstrap should be an integer. Politis and White suggest 
%       rounding the output UP to the nearest integer.
%
% Input Arguments:
%
%       data	-	[double]	numObs-by-numVars matrix, data on which we will apply stationary bootstrap or circular block bootstrap
%
%       where
%           numObs is the sample size
%           numVars is the number of variables
%
% Optional Input Arguments:
%
%       none
%
% Output Arguments:
%
%       Bstar	-	[double]    2-by-numVars vector of optimal bootstrap block lengths ([BstarSB;BstarCB])
%
% Optional Output Arguments:
%
%       none
%
% Links:
%
%       lam (subfunction), lagmat
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
%       (c) Andrew Patton, 2002
% 
%       Helpful suggestions for this code were received from Dimitris Politis and Kevin Sheppard.
%       Revisions:	To include CB - Date: 13/01/2003
%                   Speed issues (suggestions from Kevin Sheppard) - Date: 26/8/2003 
%                   Major revision: estimate of optimal block length for the stationary bootstrap corrected using 
%                       recent paper by Nordman (Annals of Statistics). - Date: 16/12/2007
%                   Fixed small bug in choice of “KN”, (thanks to Jeff Racine for spotting this) - Date: 09/12/2008
%                   Cleanup the file to remove unnecessary steps - Gabriel Bruneau - Date: 22/02/2012


% Input and Output arguments checking
% ___________________________________

narginchk(1,1); nargoutchk(0,1);


% Define size
% ___________

[numObs numVars] = size(data);


% Maximum lag selection
% _____________________

% KN extra lags to employ Politis' (2002) suggestion for finding largest significant m. We follow the empirical rule suggested 
% in Politis, 2002, "Adaptive Bandwidth Choice", as suggested in Remark 2.3, by setting c = 2, KN = 5
c = 2;
KN = max(5,sqrt(log10(numObs)));
mmax = ceil(sqrt(numObs)) + KN;             % Original: mmax = ceil(sqrt(numObs));

% Maximum value of Bstar to consider
Bmax = ceil(min(3*sqrt(numObs),numObs/3));	% December 2007: new idea for rule-of-thumb to put upper bound on estimated optimal 
origdata = data;                            % block length. Original: Bmax = sqrt(numObs)
Bstar_final = zeros(2,numVars);

for i = 1:numVars
	data = origdata(:,i);
   
	% FIRST STEP: finding mhat -> the largest lag for which the autocorrelation is still significant
    temp = lagmat(data,1:1:mmax);
	temp = temp(mmax+1:end,:);                                      % dropping the first mmax rows, as they're filled with zeros
	temp = corrcoef([data(mmax+1:end),temp]);
	temp = temp(2:end,1);

	temp2 = [lagmat(temp,1:1:KN)',temp(end-KN+1:end)];           % looking at vectors of autocorrels, from lag mhat to lag mhat+KN
    temp2 = temp2(:,KN+1:end);                                      % dropping the first KN-1, as the vectors have empty cells
	temp2 = (abs(temp2)<(c*sqrt(log10(numObs)/numObs)*ones(KN,mmax-KN+1)));	% checking which are less than the critical value
	temp2 = sum(temp2)';                                            % this counts the number of insignificant autocorrelations
	temp3 = [(1:1:length(temp2))',temp2];
	temp3 = temp3(temp2==KN,:);                                     % selecting all rows where ALL KN autocorrels are not significant
	if isempty(temp3)
        mhat = find(abs(temp) > (c*sqrt(log10(numObs)/numObs)),1,'last');     % this means that NO collection of KN autocorrels were all  
	else                                                            % insignificant, so pick largest significant lag
    	mhat = temp3(1,1);                                          % if more than one collection is possible, choose the smallest m
	end
	if 2*mhat > mmax
    	M = mmax;
	else
    	M = 2*mhat;
	end
    
	% Secon step: computing the inputs to the function for Bstar
	kk = (-M:1:M)';
   
	if M > 0
    	temp4 = lagmat(data,1:1:M);
        temp4 = temp4(M+1:end,:);               % dropping the first mmax rows, as they're filled with zeros
        temp4 = cov([data(M+1:end),temp4]);
        acv = temp4(:,1);                       % autocovariances
        acv2 = [-(1:1:M)',acv(2:end)];
        if size(acv2,1) > 1
            acv2 = sortrows(acv2,1);
        end
        acv = [acv2(:,2);acv];                  %#ok<AGROW> % autocovariances from -M to M
        Ghat = sum(lam(kk/M).*abs(kk).*acv);
        DSBhat = 2*(sum(lam(kk/M).*acv)^2);     % first part of DSBhat (note cos(0)=1) (New: After December 2007)
        %DSBhat = 2/pi*quadl('opt_block_length_calc',-pi,pi,[],[],kk,acv,lam(kk/M)); (Old: Before December 2007)
        %DSBhat = DSBhat + 4*sum(lam(kk/M).*acv)^2;     % first part of DSBhat (note cos(0)=1)
        DCBhat = 4/3*sum(lam(kk/M).*acv)^2;
      
        % FINAL STEP: constructing the optimal block length estimator
        BstarSB = ((2*(Ghat^2)/DSBhat)^(1/3))*(numObs^(1/3));
        if BstarSB > Bmax
            BstarSB = Bmax;
        end
        BstarCB = ((2*(Ghat^2)/DCBhat)^(1/3))*(numObs^(1/3));
        if BstarCB > Bmax
        	BstarCB = Bmax;
        end      
        Bstar = [BstarSB;BstarCB];
	else
        % Final step: constructing the optimal block length estimator
        Bstar = [1;1];
	end % if {M > 0}
	Bstar_final(:,i) = Bstar;
end % for {numVars}
Bstar = Bstar_final;


end % function {bootoptblock}


%--------------------------------------------------------------------------
function lam = lam(kk)	% Subfunction
% LAM:	Helper function, calculates the flattop kernel weights

lam = (abs(kk)>=0).*(abs(kk)<0.5)+2*(1-abs(kk)).*(abs(kk)>=0.5).*(abs(kk)<=1);


end % subfunction {lam}

