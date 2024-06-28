function f = exact(t,x,m)

alpha = 1/(m-1+2);
beta = alpha;
kappa = beta*(m-1)/(2*m);
k = 1;



f = (1/t^(alpha))*max(0,(k-kappa*(x.^2/(t^(2*beta))))).^(1/(m-1));

