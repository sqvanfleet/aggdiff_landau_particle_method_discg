function f = exact_2d(t,vx,vy)

d = 2;
K = 1-exp(-t/8)/2;
f = 1/(2*pi*K)^(d/2)*exp(-(vx.^2+vy.^2)/(2*K)).*(((d+2)*K-d)/(2*K)+(1-K)/(2*K^2)*(vx.^2+vy.^2));