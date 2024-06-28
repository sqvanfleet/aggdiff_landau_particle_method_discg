function f = exact(t,x)

f = (2*pi*(1-exp(-2*t)))^(-1/2)*exp(-x.^2./(2*(1-exp(-2*t))));

