%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model with a
%%%% stochastic infection rate beta
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
% note: set Refine switch to avoid interpolation
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];%optimset('MaxFunEvals',10000,'MaxIter',7000,'TolFun',...
              %                                       1e-30,'TolX',1e-30);
options3 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formatOut = 'dd-mm-yy';
DeathOld = Death;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extend = time(zz);
dates = t_span(1)+30:1:t_span(end)-extend;
Err1 = zeros(length(dates)-1,5);
Err2 = Err1;
coefMR = zeros(length(dates)-1,3);
coefMean = zeros(length(dates)-1,1);


parfor vv = 2:length(dates)

params = paramsOld;
yinit = yinitOld;
NP=1;

len0 = length(t_span(1):dates(vv))+1;
t_actualB = t_actual(1:len0);
t_extend = [t_actualB,t_actualB(end)+1:t_actualB(end)+extend];
t_spanA = t_span(1):dates(vv);
t_spanB = [t_spanA(1:end),t_spanA(end)+1:t_spanA(end)+extend];
len1 = length(t_spanB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extension of time-dependent rates of infection and death.
len = 5; s = 0;
Death = DeathOld;
Death = [Death(1:len0-s),ones(1,extend+s)*mean(Death(len0-len-s:len0-s))];
params.factorDeath = @(t)interp1(t_extend,Death,t);

%%% We shall test different extensions for the transmission parameter BETA.
%%% MEAN of last "len" days before "s" days
len = 5; s = 0;
BETA_MEAN = [BETA(1:len0-s);ones(extend+s,1)*mean(BETA(len0-len-s:len0-s))];
coefMean(vv-1) = mean(BETA(len0-len-s:len0-s));


%%% LINEAR REGRESSION of the last "len" days before "s" days
len = 10; s = 0;
t = t_actual(len0-len-s:len0-s);
OF = @(x)(BETA(len0-len-s:len0-s)-(x(1)*t + x(2)));
x = lsqnonlin(OF,[1,1]);
t =t_actual(len0-s+1:len0-s+extend)';
BETA_LR = [BETA(1:len0-s);x(1)*t + x(2)];

%%% FOURIER SERIES fitting the last "len" days before "s" days
len = 30; s = 0;
OF = @(coef)(BETA(len0-len-s:len0-s)'-FourierSeries(t_actual(len0-len-s:len0-s),coef,t_actual(len0-s)-t_actual(len0-len-s)));
coef0 = ones(1,7);
coef = lsqnonlin(OF,coef0);
t = t_actual(len0-s+1:len0-s+extend);
BETA_Fourier = [BETA(1:len0-s);FourierSeries(t,coef,t_actual(len0-s)-t_actual(len0-len-s))'];

%%% FOURIER SERIES fitting the last "len" days before "s" days with
%%% Gaussian noise
NSamples = 5000;
delta = BETA(len0-len-s:len0-s)'-FourierSeries(t_actual(len0-len-s:len0-s),coef,t_actual(len0-s)-t_actual(len0-len-s));
NOISE = std(delta)*randn(length(t),NSamples);
BETA_Fourier2 = BETA_Fourier*ones(1,NSamples);
BETA_Fourier2(len0-s+1:len0-s+extend,:) = max(0,BETA_Fourier2(len0-s+1:len0-s+extend,:)+NOISE);

%%% MEAN REVERTING Beta
dt =1/365;
ObjFun = @(coef)ObjFunHeston2(BETA(len0-len-s:len0-s),coef,t_actual(len0-len-s:len0-s),dt,BETA(len0-len-s));
coef0 = [0.1,BETA(len0-len-s),0.01];
% LB = [0,1E-3,0];
LB = [0,0,0];
UB = [10,10,10]; 
% UB = [100,10,10]; 
coef = lsqnonlin(ObjFun,coef0,LB,UB);
coefMR(vv-1,:)=coef; % [kappa,theta,xi]

beta0 = BETA(len-s+1);
kappa = coef(1);
theta = coef(2);
xi = coef(3);

NSamples = 5000;
BETA_MR = zeros(len1+1,NSamples);
BETA_MR(1:len0-s+1,:) = BETA(1:len0-s+1)*ones(1,NSamples);
beta0 = max(0,beta0)*ones(NSamples,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluating the predictions

%%% MEAN
yb = zeros(len1+1,length(yinit));
yb(1,:) = yinit;
yinit2 = yinit;
for jj = 1:len1
params.beta= BETA_MEAN(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_MEAN = sigma*yb(:,2)*N;

%%% LINEAR REGRESSION
yinit2 = yinit;
for jj = 1:len1
params.beta= BETA_LR(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_LR = sigma*yb(:,2)*N;

%%% FOURIER SERIES
yinit2 = yinit;
for jj = 1:len1
params.beta= BETA_Fourier(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_Fourier = sigma*yb(:,2)*N;

%%% FOURIER SERIES WITH NOISE

yinit2 = zeros(1,5*NSamples);
yinit2(1:NSamples) = yinit(1)*ones(1,NSamples);
yinit2(NSamples+1:2*NSamples) = yinit(2)*ones(1,NSamples);
yinit2(2*NSamples+1:3*NSamples) = yinit(3)*ones(1,NSamples);
yinit2(3*NSamples+1:4*NSamples) = yinit(4)*ones(1,NSamples);
yinit2(4*NSamples+1:5*NSamples) = yinit(5)*ones(1,NSamples);

yb = zeros(len1+1,length(yinit2));
yb(1,:) = yinit2;

for jj = 1:len1
beta = BETA_Fourier2(jj+1,:)';
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2Heston(t,y,params,beta,NSamples),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_Fourier2 = median(sigma*yb(:,NSamples+1:2*NSamples)'*N)';

%%% MEAN REVERTING
yinit2 = zeros(1,5*NSamples);
yinit2(1:NSamples) = yinit(1)*ones(1,NSamples);
yinit2(NSamples+1:2*NSamples) = yinit(2)*ones(1,NSamples);
yinit2(2*NSamples+1:3*NSamples) = yinit(3)*ones(1,NSamples);
yinit2(3*NSamples+1:4*NSamples) = yinit(4)*ones(1,NSamples);
yinit2(4*NSamples+1:5*NSamples) = yinit(5)*ones(1,NSamples);

yb = zeros(length(t_actual),length(yinit2));
yb(1,:) = yinit2;

for jj = 1:len1
tspan = t_extend(jj:jj+1);
if jj>=len0-s+1
beta = max(0,beta0 + kappa*(theta-beta0)*dt + xi*sqrt(beta0*dt).*randn(NSamples,1));
BETA_MR(jj+1,:) = beta;
else
beta = BETA_MR(jj+1,:)';    
end
[~,y]=ode45(@(t,y)seir_death_age_beta_b2Heston(t,y,params,beta,NSamples),tspan,yinit2,options);
beta0 = beta;
yinit2 = y(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_MR = median(sigma*yb(:,NSamples+1:2*NSamples)'*N)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 1;
ts = t_span(len0-s-1:len0-s+extend-1);

% %%% Evaluating Errors:
D = cumsum(data2(len0-s-1:len0-s+extend-1,1));
P1 = cumsum(NewInfections_MEAN(len0-s:len0-s+extend));
P2 = cumsum(NewInfections_LR(len0-s:len0-s+extend));
P3 = cumsum(NewInfections_Fourier(len0-s:len0-s+extend));
P4 = cumsum(NewInfections_Fourier2(len0-s:len0-s+extend));
P5 = cumsum(NewInfections_MR(len0-s:len0-s+extend));
Error = zeros(1,5);
Error(1) = abs(P1(end)-D(end))./D(end);
Error(2) = abs(P2(end)-D(end))./D(end);
Error(3) = abs(P3(end)-D(end))./D(end);
Error(4) = abs(P4(end)-D(end))./D(end);
Error(5) = abs(P5(end)-D(end))./D(end);
Err1(vv-1,:) = Error;

s=0;
D = (data2(len0-s-1:len0-s+extend-1,1));
P1 = (NewInfections_MEAN(len0-s:len0-s+extend));
P2 = (NewInfections_LR(len0-s:len0-s+extend));
P3 = (NewInfections_Fourier(len0-s:len0-s+extend));
P4 = (NewInfections_Fourier2(len0-s:len0-s+extend));
P5 = (NewInfections_MR(len0-s:len0-s+extend));

Error = zeros(1,5);
Error(1) = norm(P1-D)./norm(D);
Error(2) = norm(P2-D)./norm(D);
Error(3) = norm(P3-D)./norm(D);
Error(4) = norm(P4-D)./norm(D);
Error(5) = norm(P5-D)./norm(D);
Err2(vv-1,:) = Error;
end
