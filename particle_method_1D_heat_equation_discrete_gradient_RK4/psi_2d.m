function y = psi_2d(vx,vy,eps)

% Returns a matrix of size vx and vy with values of psi evaluated at
% (vx,vy).

d = 2;
y = 1/(2*pi*eps)^(d/2)*exp(-(vx.^2+vy.^2)/(2*eps));