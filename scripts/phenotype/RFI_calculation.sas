dm'log; clear; output; clear;odsresults;clear;';
proc datasets library=work kill memtype=data;quit;
%Let path = C:\Users\;

%let XL_in = ;
libname HFM "&path";

dm'log; clear; output; clear; odsresults; clear;';
PROC IMPORT OUT= RFI_calculation
    DATAFILE= "&path.RumenMicrobiomeGenomics_RumenNH3-N_SAS.csv" 
    DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;
proc print data= RFI_calculation (obs=5);run;
proc contents data=RFI_calculation;run;


dm'log; clear; output; clear;odsresults;clear;';
ods trace on;
ods select Tests3 Fitstatistics Contrasts LSMeans Diffs Estimates;
ods excel file = "&path.RFI_output.xlsx";
proc mixed data=RFI_calculation method=type3 covtest ratio itdetails maxiter=500;
class ExpN Parity;
model DMI = Parity NESec MBW BEC  /ddfm=kr s residual outp=residuals;
random ExpN / Solution;
lsmeans Parity / pdiff adjust=Tukey adjdfe=Row;
ods output rs=resid output=pred;
run;
ods select all;
ods trace off;
ods excel close;


ods excel file = "&path.RFI_output.xlsx";
proc print data=pred;
run;
ods excel close;
