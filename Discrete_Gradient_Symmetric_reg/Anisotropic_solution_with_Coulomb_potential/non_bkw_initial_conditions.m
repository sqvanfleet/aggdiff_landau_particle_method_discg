function f = non_bkw_initial_conditions(vx,vy)

u_1 = [-2,1]; u_2 = [0,-1];

f = (1/4)*pi*(exp(-((vx-u_1(1)).^2+(vy-u_1(2)).^2)/2)+exp(-((vx-u_2(1)).^2+(vy-u_2(2)).^2)/2));

% function f = f_initial_2D_Coulumb(vx,vy)
% 
% f = (1/(4*pi))*(exp(-((vx+2).^2+(vy-1).^2)/2)+exp(-(vx.^2+(vy+1).^2)/2));