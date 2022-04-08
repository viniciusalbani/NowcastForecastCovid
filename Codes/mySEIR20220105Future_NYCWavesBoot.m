clear all; clc; close all; format long e;
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
% load('data_NYC_20220103Boot.mat')
load('data_NYC_20220328Boot.mat')
formatOut = 'dd-mmm-yy';
DeathOld = Death;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extend = 45;
dates = [datetime(2020,04,05),datetime(2020,10,18),datetime(2021,07,10),datetime(2021,11,01)];
Err1 = zeros(length(dates)-1,5);
Err2 = Err1;
H = [100,100,700,400];
NSamples0 = NSamples;

for vv = length(dates)
    if length(dates(end):t_span(end))<extend
        extend = length(dates(end):t_span(end))-1;
    end
len0 = length(t_span(1):dates(vv))+1;
t_actualB = t_actual(1:len0);
t_extend = [t_actualB,t_actualB(end)+1:t_actualB(end)+extend];
t_spanA = t_span(1):dates(vv);
t_spanB = [t_spanA(1:end),t_spanA(end)+1:t_spanA(end)+extend];
len1 = length(t_spanB);
    
NewInfections_MEAN = zeros(len1+1,NSamples0);
NewInfections_LR = zeros(len1+1,NSamples0);
NewInfections_Fourier = zeros(len1+1,NSamples0);
NewInfections_Fourier2 = zeros(len1+1,NSamples0);
NewInfections_MR = zeros(len1+1,NSamples0);
parfor zz = 1:NSamples0
params = paramsOld;
yinit = yinitOld;
NP=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extension of time-dependent rates of infection and death.
len = 5; s = 0;
Death = DeathOld;
Death = [Death(1:len0-s),ones(1,extend+s)*mean(Death(len0-len-s:len0-s))];
params.factorDeath = @(t)interp1(t_extend,Death,t);

%%% We shall test different extensions for the transmission parameter BETA.
%%% MEAN of last "len" days before "s" days
len = 5; s = 0;
BETA_MEAN = [BETABoot(zz,1:len0-s)';ones(extend+s,1)*mean(BETABoot(zz,len0-len-s:len0-s))'];

%%% LINEAR REGRESSION of the last "len" days before "s" days
len = 10; s = 0;
t = t_actual(len0-len-s:len0-s);
OF = @(x)(BETABoot(zz,len0-len-s:len0-s)-(x(1)*t + x(2)));
x = lsqnonlin(OF,[1,1]);
t =t_actual(len0-s+1:len0-s+extend)';
BETA_LR = [BETABoot(zz,1:len0-s)';x(1)*t + x(2)];

%%% FOURIER SERIES fitting the last "len" days before "s" days
len = 30; s = 0;
OF = @(coef)(BETABoot(zz,len0-len-s:len0-s)-FourierSeries(t_actual(len0-len-s:len0-s),coef,t_actual(len0-s)-t_actual(len0-len-s)));
coef0 = ones(1,7);
coef = lsqnonlin(OF,coef0);
t = t_actual(len0-s+1:len0-s+extend);
BETA_Fourier = [BETABoot(zz,1:len0-s)';FourierSeries(t,coef,t_actual(len0-s)-t_actual(len0-len-s))'];

%%% FOURIER SERIES fitting the last "len" days before "s" days with
%%% Gaussian noise
NSamplesB = 5000;
delta = BETABoot(zz,len0-len-s:len0-s)-FourierSeries(t_actual(len0-len-s:len0-s),coef,t_actual(len0-s)-t_actual(len0-len-s));
NOISE = std(delta)*randn(length(t),NSamplesB);
BETA_Fourier2 = BETA_Fourier*ones(1,NSamplesB);
BETA_Fourier2(len0-s+1:len0-s+extend,:) = max(0,BETA_Fourier2(len0-s+1:len0-s+extend,:)+NOISE);

%%% MEAN REVERTING Beta
% len = 30; s = 0;
dt =1/365;
ObjFun = @(coef)ObjFunHeston2(BETABoot(zz,len0-len-s:len0-s)',coef,t_actual(len0-len-s:len0-s),dt,BETABoot(zz,len0-len-s));
coef0 = [0.1,BETABoot(zz,len0-len-s),0.01];
LB = [0,0,0];
UB = [10,10,10];
coef = lsqnonlin(ObjFun,coef0,LB,UB);

beta0 = BETABoot(zz,len-s+1);
kappa = coef(1);
theta = coef(2);
xi = coef(3);

NSamplesB = 5000;
BETA_MR = zeros(len1+1,NSamplesB);
BETA_MR(1:len0-s+1,:) = BETABoot(zz,1:len0-s+1)'*ones(1,NSamplesB);
beta0 = max(0,beta0)*ones(NSamplesB,1);

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
NewInfections_MEAN(:,zz) = sigma*yb(:,2)*N;

%%% LINEAR REGRESSION
yinit2 = yinit;
for jj = 1:len1
params.beta= BETA_LR(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_LR(:,zz) = sigma*yb(:,2)*N;

%%% FOURIER SERIES
yinit2 = yinit;
for jj = 1:len1
params.beta= BETA_Fourier(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections_Fourier(:,zz) = sigma*yb(:,2)*N;

%%% FOURIER SERIES WITH NOISE

yinit2 = zeros(1,5*NSamplesB);
yinit2(1:NSamplesB) = yinit(1)*ones(1,NSamplesB);
yinit2(NSamplesB+1:2*NSamplesB) = yinit(2)*ones(1,NSamplesB);
yinit2(2*NSamplesB+1:3*NSamplesB) = yinit(3)*ones(1,NSamplesB);
yinit2(3*NSamplesB+1:4*NSamplesB) = yinit(4)*ones(1,NSamplesB);
yinit2(4*NSamplesB+1:5*NSamplesB) = yinit(5)*ones(1,NSamplesB);

yb = zeros(len1+1,length(yinit2));
yb(1,:) = yinit2;

for jj = 1:len1
beta = BETA_Fourier2(jj+1,:)';
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2Heston(t,y,params,beta,NSamplesB),t_extend(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections = sigma*yb(:,NSamplesB+1:2*NSamplesB)*N;

% aux = sort(NewInfections');
% aux2 = round(0.15*NSamples);
% aux = aux(aux2+1:end-aux2,:);
% CI90FourierWN = [min(aux);max(aux)];
NewInfections_Fourier2(:,zz) = median(NewInfections');

%%% MEAN REVERTING
yinit2 = zeros(1,5*NSamplesB);
yinit2(1:NSamplesB) = yinit(1)*ones(1,NSamplesB);
yinit2(NSamplesB+1:2*NSamplesB) = yinit(2)*ones(1,NSamplesB);
yinit2(2*NSamplesB+1:3*NSamplesB) = yinit(3)*ones(1,NSamplesB);
yinit2(3*NSamplesB+1:4*NSamplesB) = yinit(4)*ones(1,NSamplesB);
yinit2(4*NSamplesB+1:5*NSamplesB) = yinit(5)*ones(1,NSamplesB);

yb = zeros(len1+1,length(yinit2));
yb(1,:) = yinit2;

for jj = 1:len1
tspan = t_extend(jj:jj+1);
if jj>=len0-s+1
beta = max(0,beta0 + kappa*(theta-beta0)*dt + xi*sqrt(beta0*dt).*randn(NSamplesB,1));
BETA_MR(jj+1,:) = beta;
else
beta = BETA_MR(jj+1,:)';    
end
[~,y]=ode45(@(t,y)seir_death_age_beta_b2Heston(t,y,params,beta,NSamplesB),tspan,yinit2,options);
beta0 = beta;
yinit2 = y(end,:);
yb(jj+1,:) = yinit2;
end 
NewInfections = sigma*yb(:,NSamplesB+1:2*NSamplesB)*N;

aux = sort(NewInfections');
aux2 = round(0.25*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI90MR = [min(aux);max(aux)];
NewInfections_MR(:,zz)= median(NewInfections');
CI901MR(:,zz) = CI90MR(1,:)';
CI902MR(:,zz) = CI90MR(2,:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 1;
tt = [10,15,30,45,60];
ts = t_span(len0-s:len0-s+extend);

%%% Evaluating Errors:
D = cumsum(data2(len0-s-1:len0-s+extend-1,1));

aux = sort(NewInfections_MEAN(len0-s:len0-s+extend,:)');
aux2 = round(0.05*NSamples0);
aux = aux(aux2+1:end-aux2,:);
CI90Mean = cumsum([min(aux);max(aux)]')';
P1 = cumsum(median(NewInfections_MEAN(len0-s:len0-s+extend,:)'));

aux = sort(NewInfections_LR(len0-s:len0-s+extend,:)');
aux2 = round(0.05*NSamples0);
aux = aux(aux2+1:end-aux2,:);
CI90LR = cumsum([min(aux);max(aux)]')';
P2 = cumsum(median(NewInfections_LR(len0-s:len0-s+extend,:)'));

aux = sort(NewInfections_Fourier(len0-s:len0-s+extend,:)');
aux2 = round(0.05*NSamples0);
aux = aux(aux2+1:end-aux2,:);
CI90Fourier = cumsum([min(aux);max(aux)]')';
P3 = cumsum(median(NewInfections_Fourier(len0-s:len0-s+extend,:)'));

aux = sort(NewInfections_Fourier2(len0-s:len0-s+extend,:)');
aux2 = round(0.05*NSamples0);
aux = aux(aux2+1:end-aux2,:);
CI90Fourier2 = cumsum([min(aux);max(aux)]')';
P4 = cumsum(median(NewInfections_Fourier2(len0-s:len0-s+extend,:)'));

aux = sort(NewInfections_MR(len0-s:len0-s+extend,:)');
aux2 = round(0.05*NSamples0);
aux = aux(aux2+1:end-aux2,:);
CI90MR = cumsum([min(aux);max(aux)]')';
P5 = cumsum(median(NewInfections_MR(len0-s:len0-s+extend,:)'));

figure
hold on
box on
title('Average')
h1=area(ts,CI90Mean(2,:),'linestyle',':','FaceColor','r','FaceAlpha',0.5);
h2=area(ts,CI90Mean(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(ts,D,'FaceColor',[0 0.75 0.75],'EdgeColor','k','FaceAlpha',0.5)
plot(ts,P1,'r','LineWidth',2')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Observed','Predicted','Location','NorthWest')
xlim([ts(1)-1,ts(end)])
% xticks([0,tt])
ylabel('Accum. Infections')
% xlabel('Number of Days Ahead')
set(gcf,'Position',H)
set(gca,'FontSize',22,'FontName','Arial')
hold off
saveas(gcf,['NYC_Mean',num2str(vv),'.fig']);
print('-dpng',['NYC_Mean',num2str(vv)]);

figure
hold on
box on
title('Linear Regression')
h1=area(ts,CI90LR(2,:),'linestyle',':','FaceColor','r','FaceAlpha',0.5);%[51,236,255]/255);
h2=area(ts,CI90LR(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(ts,D,'FaceColor',[0 0.75 0.75],'EdgeColor','k','FaceAlpha',0.5)
plot(ts,P2,'r','LineWidth',2')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Observed','Predicted','Location','NorthWest')
xlim([ts(1)-1,ts(end)])
ylabel('Accum. Infections')
set(gcf,'Position',H)
set(gca,'FontSize',22,'FontName','Arial')
hold off
saveas(gcf,['NYC_LR',num2str(vv),'.fig']);
print('-dpng',['NYC_LR',num2str(vv)]);


figure
hold on
box on
title('Fourier')
h1=area(ts,CI90Fourier(2,:),'linestyle',':','FaceColor','r','FaceAlpha',0.5);%[51,236,255]/255);
h2=area(ts,CI90Fourier(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(ts,D,'FaceColor',[0 0.75 0.75],'EdgeColor','k','FaceAlpha',0.5)
plot(ts,P3,'r','LineWidth',2')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Observed','Predicted','Location','NorthWest')
xlim([ts(1)-1,ts(end)])
ylabel('Accum. Infections')
set(gcf,'Position',H)
set(gca,'FontSize',22,'FontName','Arial')
hold off
saveas(gcf,['NYC_Fourier',num2str(vv),'.fig']);
print('-dpng',['NYC_Fourier',num2str(vv)]);

figure
hold on
box on
title('Fourier with Noise')
h1=area(ts,CI90Fourier2(2,:),'linestyle',':','FaceColor','r','FaceAlpha',0.5);%[51,236,255]/255);
h2=area(ts,CI90Fourier2(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(ts,D,'FaceColor',[0 0.75 0.75],'EdgeColor','k','FaceAlpha',0.5)
plot(ts,P4,'r','LineWidth',2')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Observed','Predicted','Location','NorthWest')
xlim([ts(1)-1,ts(end)])
ylabel('Accum. Infections')
set(gcf,'Position',H)
set(gca,'FontSize',22,'FontName','Arial')
hold off
saveas(gcf,['NYC_Fourier2',num2str(vv),'.fig']);
print('-dpng',['NYC_Fourier2',num2str(vv)]);

figure
hold on
box on
title('Mean Reverting')
h1=area(ts,CI90MR(2,:),'linestyle',':','FaceColor','r','FaceAlpha',0.5);%[51,236,255]/255);
h2=area(ts,CI90MR(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(ts,D,'FaceColor',[0 0.75 0.75],'EdgeColor','k','FaceAlpha',0.5)
plot(ts,P5,'r','LineWidth',2')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Observed','Predicted','Location','NorthWest')
xlim([ts(1)-1,ts(end)])
ylabel('Accum. Infections')
set(gcf,'Position',H)
set(gca,'FontSize',22,'FontName','Arial')
hold off
saveas(gcf,['NYC_MR',num2str(vv),'.fig']);
print('-dpng',['NYC_MR',num2str(vv)]);
end
