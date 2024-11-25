function y = psi_1d(x,eps)

d = 1;
y = 1/(2*pi*eps)^(d/2)*exp(-x.^2/(2*eps));