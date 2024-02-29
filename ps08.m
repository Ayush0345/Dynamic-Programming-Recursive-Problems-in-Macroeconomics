% upload the data here
% use readnatrix command
% at the end of this block of code you should have 3 variables: y, c and
% inv

y =  readmatrix("GDPC1.xls") ;
c = readmatrix("PCECC96.xls") ;
inv = readmatrix("GPDIC1.xls");
weekly = readmatrix("PRS85006023.xls");
emp = readmatrix("CE16OV.xls");

% =========================================================================
% Connstruct N = weekly x emp

N = weekly*emp(307,1).';

% =========================================================================
% Data manipulation:  take logs and apply the HP filter

% Apply the HP filter with smoothing parameter 1600 to the logged variables
% For a generic variable x, call log x by xlog
% keep the cycle component and trend component
% call the trend component of xlog  xlogtrend; 
% call the cyclical component of xlog  xlogcycle 

ylog = log(y);
clog = log(c);
invlog = log(inv);
wlog = log(weekly);
elog = log(emp);
nlog = log(N);
    
% many variables here!

lambda = 1600;

% y
[ytrend, ycyclical] = hpfilter(ylog, lambda);

% c
[ctrend, ccyclical] = hpfilter(clog, lambda);

% inv
[itrend, icyclical] = hpfilter(invlog, lambda);

% N
[ntrend, ncyclical] = hpfilter(nlog, lambda);

% ========================================================================
% plot the data and the smoothed (trend component) of ylog; check command
% hold on


 .... Add code here.... 
plot(ylog, ytrend, "Red");
hold on
axis padded 
title("Log of Real GDP v/s Detrended GDP")
xlabel("Log of Real GDP")
ylabel("Detrended GDP")

% ========================================================================
% Find standard deviations
% use std command to find the standard deviation of the cycle components 
% of all variables 


 .... Add code here.... 
 Ystd = std(ycyclical);
 Cstd = std(ccyclical);
 Istd = std(icyclical);
 Nstd = std(ncyclical);


% Add code to display results; check command disp or command fprintf

disp(Ystd)
disp(Cstd)
disp(Istd)
disp(Nstd)
% ========================================================================
% find relative standard deviation 

% Calculating mean to measure relative volatility

mean_y = mean(ycyclical);
mean_c = mean(ccyclical);
mean_i = mean(icyclical);
mean_n = mean(ncyclical);

% .... Add code here.... 

% Relative Std = Std/Mean

Ystd2 = Ystd/mean_y;
Cstd2 = Cstd/mean_c;
Istd2 = Istd/mean_i;
Nstd2 = Nstd/mean_n;

% .... Add code to display results here.... 

disp(Ystd2)
disp(Cstd2)
disp(Istd2)
disp(Nstd2)

% ========================================================================
% Look at correlations
% use the corrcoef function to find the correlation of xlogcycle with ylogcycle
R = corrcoef(ycyclical, ccyclical);
P = corrcoef(ycyclical, icyclical);
K = corrcoef(ycyclical, ncyclical);

disp(R)
disp(P)
disp(K)

% Part (i) - The standard deviation of y, c, i, and n are slightly less 
% than the values in the business cycle statistics from the paper. This is
% because the data we are analysing is till quarter 3 of 2023 and the data
% analysed in the paper is much older (2000). Therefore, the relative
% volatility of c, i, and n is lesser than the values in the paper because
% the volatility has smoothened over time with less dispersion of risk
% (variance) in the data in the next 23 years. A major shock in the
% business cycle came during the 2008 global financial crisis which has
% been estimated from the data cumulatively.

% ========================================================================


