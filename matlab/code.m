addpath 'C:\Users\andreas.timoudas\Python\MasterThesis\matlab\MATLAB_TVARX_Toolbox\'
addpath 'C:\Users\andreas.timoudas\Python\MasterThesis\matlab\MATLAB_TVARX_Toolbox\utilities\'
import MATLAB_TRARX_Toolbox.*

table = readtable('data.csv');
y_data = table(:,[3 4 5 6]);
y = table2array(y_data);

options_ic = struct('infoc', 'AIC');
options_tvarx = struct('optimcrit', 'loglik');
th = logical([1 0 0 0]);
vnames = char('SFSI', 'BNP', 'CPI', 'DINT');

%[icVar, infocrit, tresh] = tvarxic(y,10,1,th, options_ic, options_tvarx);

result = tvarx(y, 3, 1, th, 2, struct('prtr', true, 'optimcrit', 'logdet', 'const', true));

irf_options = struct('sqrttype', 'chol', 'vnames', vnames, 'graph', true); 
%[IRF,bootIRFLowerCI,bootIRFUpperCI,bootIRFMedCI,bootIRFAll,impact] = tvarxlirf(result, irf_options);



