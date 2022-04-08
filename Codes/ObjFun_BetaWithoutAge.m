function f = ObjFun_BetaWithoutAge(t_actual,params,data,options,priors,yinit,unknowns)
Number = params.NumberOfAgeClasses;
params.beta = unknowns(1);
tspan = [t_actual(1),t_actual(end)];

N = params.N;


sigma = params.sigma;

[~,y]=ode45(@(t,y)seir_death_age_beta_b2(t,y, params),tspan,yinit,options);
NewInfections = sigma*y(end,Number+1:2*Number)*N;

% %%% Gaussian Misfit or Likelihood
% f = NewInfections-data(:,1);

%%% log-Poisson Misfit or Likelihood with Stirling Formula 
Stirling = 0.5*log(2*pi*data(:,1)) + data(:,1).*log(data(:,1)) - data(:,1);
f = data(:,1).*log(NewInfections) - NewInfections - Stirling;

%%% Including Priors
f = [f;1E-3*(unknowns-priors)'];
f(isnan(f)~=0)=zeros;