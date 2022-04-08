function f = FourierSeries(x,coef,L)
n = (length(coef)-1)/2;
f = coef(1)*ones(size(x));
aux = 2*pi/L;
for jj = 1:n
    f = f+coef(2*jj)*cos(aux*jj*x) + coef(2*jj+1)*sin(aux*jj*x);
end