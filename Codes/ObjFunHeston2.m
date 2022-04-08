function f = ObjFunHeston2(data,coef,t,dt,beta0)
kappa = coef(1);
theta = coef(2);
xi = coef(3);

kappa2 = kappa^2;
theta2 = theta^2;
xi2 = xi^2;
data2 = data.^2;
dt2 = dt^2;

MeanBeta = beta0;
MeanBetaSqr = beta0^2;
f = zeros(size(t));
f(1) = sqrt(MeanBetaSqr -2*data(1)*MeanBeta + data2(1));
for jj = 1:length(t)-1
MeanBeta = MeanBeta + kappa*(theta-MeanBeta)*dt;
MeanBetaSqr = MeanBetaSqr + kappa2*(theta2 - 2*theta*MeanBeta + MeanBetaSqr)*dt2 + xi2*MeanBeta*dt + 2*kappa*(theta*MeanBeta - MeanBetaSqr)*dt;
f(jj) = sqrt(abs(MeanBetaSqr -2*data(1)*MeanBeta + data2(jj)));
end
f = [f,1E-3*coef];