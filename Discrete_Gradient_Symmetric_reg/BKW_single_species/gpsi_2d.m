function [y1,y2] = gpsi_2d(vx,vy,eps)

y1 = psi_2d(vx,vy,eps).*(-vx/eps);
y2 = psi_2d(vx,vy,eps).*(-vy/eps);