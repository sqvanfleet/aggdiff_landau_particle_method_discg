clc
clear all
close all


% 2D Maxwell molecule 
gamma = -3;
C_gamma = 1/16;

% computational domain
Vmax = 10;

n_list = [40,45,50,55,60];

for alpha = 1:length(n_list)

n = n_list(alpha);
% total number of particles
Np = n^2;

% initial mesh
dv = 2*Vmax/n;
v = (-Vmax+dv/2):dv:(Vmax-dv/2);
[vx,vy] = ndgrid(v);
Vx = vx(:);
Vy = vy(:);
 
% initial weight
t0 = 0; 
f0 = non_bkw_initial_conditions(vx,vy);
w = dv^2*f0;
W = w(:);

% reconstruction mesh
Nr = n;
dvr = 2*Vmax/Nr;
vr = (-Vmax+dvr/2):dvr:(Vmax-dvr/2); 
[vrx,vry] = ndgrid(vr);
Vrx = vrx(:);
Vry = vry(:);

% reconstructed solution
f = zeros(Nr,Nr);
% choosing epsilon
epsilon = 4*(0.4*(dv)^0.99)^2;
% parfor for i
parfor i = 1:Nr
    for j = 1:Nr
        f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
    end
end

% initial comparison
f_initial = non_bkw_initial_conditions(vrx,vry);
% figure(1)
% plot(vr,f_initial(:,Nr/2),'b')
% hold on
% plot(vr,f(:,Nr/2),'r')
% hold off
% title('Exact solution (blue) and par%ticle method (red)')

dt = 0.1/2;
tmax = 20;

%Nt = 10;
Nt = round((tmax-t0)/dt);
error_list = zeros(Nt,10);

tic

% max iterations and tolerence for fixed point iteration
max_iter = 300;
tol = 1e-15;
Fixed_Point_list = zeros(Nt,2);

for nt = 1:Nt
    time = t0+dt*nt;
    %Initial Guess forward Euler step
    [U_x,U_y,gF_x,gF_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx,Vy,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
    Vx_new = Vx + dt*U_x;
    Vy_new = Vy + dt*U_y;
    %Fixed point iteration to solve non-linear system that results from
    %implicit scheme

    for i = 1:max_iter
        Vx_old = Vx_new;
        Vy_old = Vy_new;
        [U_x,U_y,gF_x,gF_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx_new,Vy_new,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
        Vx_new = Vx + dt*U_x;
        Vy_new = Vy + dt*U_y;
        absres = sqrt(sum((Vx_new - Vx_old).^2 + (Vy_new - Vy_old).^2));
        relres = absres/(sqrt(sum(Vx_new.^2 + Vy_new.^2)));
        %disp(['abs res = ',num2str(absres)])
        disp(['rel res = ',num2str(relres)])
        if relres < tol
            disp(['fixed point iteration took ',num2str(i),' iterations at time ', num2str(time)])
            numiter = i;
            break
        end
        if i == max_iter
            Fixed_Point_list(nt,1) = time;
            Fixed_Point_list(nt,2) = res;
            disp(['maximum number of iterations reached at time ',num2str(time)])
            numiter = i;
        end
    end

    Vx = Vx_new;
    Vy = Vy_new;
    error_list(nt,10) = numiter;

    %tic
    % reconstructed solution
    % parfor for i
    parfor i = 1:Nr
        for j = 1:Nr
            f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
        end
    end
    %toc
    
 % plot

%     figure(2)
%     plot(vr,f(:,Nr/2),'r')
%     hold on 
%     plot(vr,f_initial(:,Nr/2),'b')
%     hold off
%     title('Particle method f(:,Nr/2)')
%     drawnow
    
    % store
    error_list(nt,1) = time;
    disp("current time: ");
    disp(time);
    
    % relative error

    
    % moments
    rho = sum(W);
    m1 = sum(W.*Vx);
    m2 = sum(W.*Vy);
    E = sum(W.*(Vx.^2+Vy.^2));   
    F = zeros(Nr^2,1);
    % parfor for i
    parfor i = 1:Nr^2 %compute the entropy at time t
        F(i) = sum(W.*psi_2d(Vrx(i)-Vx,Vry(i)-Vy,epsilon))
    end

    %compute the fisher information
   
    Fish = sum(W.*(gF_x.^2 + gF_y.^2)); %fisher information 
    eta = dv^2*sum(F.*log(F));
    error_list(nt,2) = rho;   
    error_list(nt,3) = m1;   
    error_list(nt,4) = m2;
    error_list(nt,5) = E;  
    error_list(nt,6) = eta;  
    error_list(nt,8) = Fish;
    error_list(nt,9) = dissipation;

end

total_time = toc;

filename = ['Data/particle_2d_n',num2str(n),'_dv',num2str(dv),...
    '_dt_',num2str(dt),'_time_',num2str(time),'.mat'];
save(filename)

end
