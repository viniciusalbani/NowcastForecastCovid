clear all; clc; close all; format long e;
DATA_BC;
mySEIR20220328Stochastic;
save('data_BC_20220328.mat');

clear all; clc; close all; format long e;
DATA_SP2;
mySEIR20220328Stochastic;
save('data_SP_20220328.mat');

clear all; clc; close all; format long e;
DATA_NYC;
mySEIR20220328Stochastic;
save('data_NYC_20220328.mat');

clear all; clc; close all; format long e;
DATA_Chicago;
mySEIR20220328Stochastic;
save('data_Chicago_20220328.mat');
 

clear all; clc; close all; format long e;
load('data_BC_20220328.mat')
time = [10,15,30,45,60,90];
for zz = 1:length(time)
mySEIR20220109Future;
save(['data_BC_20220101_',num2str(time(zz)),'dias.mat']);
end


clear all; clc; close all; format long e;
load('data_NYC_20220328.mat')
time = [10,15,30,45,60,90];
for zz = 3:length(time)
mySEIR20220109Future;
save(['data_NYC_20220101_',num2str(time(zz)),'dias.mat']);
end


clear all; clc; close all; format long e;
load('data_Chicago_20220328.mat')
time = [10,15,30,45,60,90];
for zz = 1:length(time)
mySEIR20220109Future;
save(['data_Chicago_20220101_',num2str(time(zz)),'dias.mat']);
end


clear all; clc; close all; format long e;
load('data_SP_20220328.mat');
time = [10,15,30,45,60,90];
for zz = 1%:length(time)
mySEIR20220109Future;
save(['data_SP_20220101_',num2str(time(zz)),'dias.mat']);
end

clc;
SummaryStatBC;
SummaryStatNYC;
SummaryStatChicago;
SummaryStatSP;