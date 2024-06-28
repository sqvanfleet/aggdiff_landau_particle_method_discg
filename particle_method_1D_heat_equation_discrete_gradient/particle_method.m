clc
clear all
close all 


n_list = [60,70,80,90,100]; %number of particles

for alpha = 1:length(n_list)

% computational domain
Xmax = 15;
% total number of particles
n = n_list(alpha);

% initial mesh
dx = 2*Xmax/n;
x = (-Xmax+dx/2):dx:(Xmax-dx/2);


t0 = 2; % initial time
f0 = exact(t0,x); %initial condition
w = dx*f0; %midpoint


% reconstruction mesh
Nr = n;
dxr = 2*Xmax/n; %cell length
xr = (-Xmax+dxr/2):dxr:(Xmax-dxr/2);  %cell centers

% reconstructed solution
f = zeros(1,Nr); %initialization


% choosing epsilon
epsilon = 4*(0.4*(dx)^0.99)^2;
% will this choice in epsilon work for the wave equation? It worked for the
% Landau Equation...

% parfor i
for i = 1:Nr
    f(i) = sum(w.*psi_1d(x-xr(i),epsilon));
end



% % initial comparison uncomment for plot
f_exact = exact(t0,x);
% figure(1)
% plot(xr,f_exact,'b','DisplayName','Initial Conditions')
% hold on
% plot(xr,f,'r','DisplayName','Particle Method Initialization')
% hold off
% legend
% title('Initialization')


tmax = 3;
dt = .01;
Nt = round((tmax-t0)/dt);
error_list = zeros(Nt,7);


max_iter = 300;
tol = 1e-15;
Fixed_Point_list = zeros(Nt,2);


for nt = 1:Nt

    time = t0+dt*nt;

    %initial forward Euler step
    dF = right_hand_side(w,x,x,xr,dx,epsilon,n);
    xnew = x - dt*dF;

    %Fixed point iteration to solve non-linear system that results from
    %implicit scheme

    for i = 1:max_iter
        x_old = xnew;
        gF = right_hand_side(w,x,xnew,xr,dx,epsilon,n);
        xnew = x - dt*gF;
        absres = sqrt(sum((xnew - x_old).^2));
        relres = absres/(sqrt(sum(xnew.^2)));
        disp(['rel res = ',num2str(relres)]);

        if relres < tol
            disp(['fixed point iteration took ',num2str(i),' iterations at time ', num2str(time)])
            numiter = i;
            break
        end
        if i == max_iter
            Fixed_Point_list(nt,1) = time;
            Fixed_Point_list(nt,2) = relres;
            disp(['maximum number of iterations reached at time ',num2str(time)])
            numiter = i;
        end
    end


    x = xnew;

    error_list(nt,7) = numiter;



    % reconstruct the particle method solution
    %parfor i 
    for i = 1:Nr
        f(i) = sum(w.*psi_1d(x-xr(i),epsilon));
    end

    disp(['Current time = ',num2str(time)])
    f_exact = exact(time,xr);

% %    figure(2) uncomment to see figure
%     plot(xr,f_exact,'b','DisplayName','Exact Solution');
%     hold on 
%     plot(xr,f,'r','DisplayName','Particle Method Solution');
%     hold off
%     xlabel('x')
%     ylabel('$\rho$','Interpreter','latex')
%     legend
%     drawnow


    % things to store
    error_list(nt,1) = time;
    f_error = f-f_exact;

    %relative error
    Linf_error = max(abs(f_error))/max(abs(f_exact));
    L1_error = sum(abs(f_error))/sum(abs(f_exact));
    L2_error = sqrt(sum(f_error.^2)/sum(f_exact.^2));
    error_list(nt,2) = Linf_error;
    error_list(nt,3) = L1_error;
    error_list(nt,4) = L2_error;


    rho = sum(w);
    inside = zeros(1,n);
    for i = 1:n
        inside(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
    end
    eta = dx*sum(inside.*log(inside));
    error_list(nt,5) = rho;
    error_list(nt,6) = eta;

end


filename = ['Data/particle_2d_n',num2str(n),'_dx',num2str(dx),'_epsilon',...
    num2str(epsilon),'_dt',num2str(dt),'_tmax',num2str(time),'.mat'];
save(filename)

end