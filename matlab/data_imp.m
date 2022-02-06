[data,names]=xlsread('datab.xls');

datestrings=names(2:end,1);
startdate = '1995q4';
enddate = '2021q3';
startlocationData=find(strcmp(datestrings,startdate));
endlocation=find(strcmp(datestrings,enddate));
data=data(startlocationData:endlocation);
