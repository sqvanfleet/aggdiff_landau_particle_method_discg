clc
clear all



% 2D Maxwell molecule 
gamma = 0;
C_gamma = 1/16;
n_list = [40,45,50,55,60];

for alpha = 1:length(n_list)
% computational domain
Vmax = 4;
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
f0 = exact_2d(t0,vx,vy);
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
f_exact = exact_2d(t0,vrx,vry);
% figure(1)
% plot(vr,f_exact(:,Nr/2),'b')
% % hold on
% % plot(vr,f(:,Nr/2),'r')
% title('Exact solution (blue) and par%ticle method (red)')


tmax = 5;
dt = 0.01/8;
%Nt = 2;
Nt = round((tmax-t0)/dt);
error_list = zeros(Nt,13);

tic

% max iterations and tolerence for fixed point iteration
max_iter = 400;
tol = 1e-15;
Fixed_Point_list = zeros(Nt,2);

for nt = 1:Nt
    time = t0+dt*nt;
    %Initial Guess forward Euler step
    [U_x,U_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx,Vy,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
    Vx_new = Vx + dt*U_x;
    Vy_new = Vy + dt*U_y;
    %Fixed point iteration to solve non-linear system that results from
    %implicit scheme
    tic
    for i = 1:max_iter
        Vx_old = Vx_new;
        Vy_old = Vy_new;
        [U_x,U_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx_new,Vy_new,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
        Vx_new = Vx + dt*U_x;
        Vy_new = Vy + dt*U_y;
        absres = sqrt(sum((Vx_new - Vx_old).^2 + (Vy_new - Vy_old).^2));
	    relres = absres/sqrt(sum(Vx_new.^2+Vy_new.^2));
        disp(['relres = ',num2str(relres)])
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
    itertime = toc;
    disp(num2str(itertime))
    Vx = Vx_new;
    Vy = Vy_new;
    error_list(nt,13) = numiter;
    %toc
    
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

   f_exact = exact_2d(time,vrx,vry);
%    figure(2)
%    plot(vr,f_exact(:,Nr/2),'b');
%    hold on
%    plot(vr,f(:,Nr/2),'r')
%    title('Exact solution (blue) and particle method (red)')
%    hold off
%    drawnow
    
    % store
    error_list(nt,1) = time;
    disp("current time: ");
    disp(time);
    
    f_error = f-f_exact;
    
    % relative error
    Linf_error = max(max(abs(f_error)))/max(max(abs(f_exact)));
    L1_error = sum(sum(abs(f_error)))/sum(sum(abs(f_exact)));
    L2_error = sqrt(sum(sum(f_error.^2))/sum(sum(f_exact.^2)));   
    error_list(nt,2) = Linf_error;   
    error_list(nt,3) = L1_error;   
    error_list(nt,4) = L2_error;
    
    % moments
    rho = sum(W);
    m1 = sum(W.*Vx);
    m2 = sum(W.*Vy);
    E = sum(W.*(Vx.^2+Vy.^2));   
    F = zeros(Nr^2,1);
    % parfor for i
    parfor i = 1:Nr^2
        F(i) = sum(W.*psi_2d(Vrx(i)-Vx,Vry(i)-Vy,epsilon));
    end
    %compute the fisher information 
    term1_x = zeros(Np,1);
    term1_y = zeros(Np,1);
    for i = 1:Np %compute fisher information and dissipation terms
        % at time t = t0 + nt*dt
        [A,B] = gpsi_2d(Vx(i)-Vrx,Vy(i)-Vry,epsilon);
        term1_x(i) = dv^2*sum(A.*log(F));
        term1_y(i) = dv^2*sum(B.*log(F));

    end
    Fish = sum(W.*(term1_x.^2 + term1_y.^2)); %fisher information 
    eta = dv^2*sum(F.*log(F));
    eta1 = dv^2*sum(F.*(log(F)+log(2*pi)+(Vrx.^2+Vry.^2)/2));
    error_list(nt,5) = rho;   
    error_list(nt,6) = m1;   
    error_list(nt,7) = m2;
    error_list(nt,8) = E;  
    error_list(nt,9) = eta;  
    error_list(nt,10) = eta1;
    error_list(nt,11) = Fish;
    error_list(nt,12) = dissipation;
end

r = Fixed_Point_list(:,2)>tol;
if max(r) == 1
    disp('Fixed point iteration achieved the maximum number of iterations at time ')
    disp(Fixed_Point_list(r,:))
end

if max(r) == 0
    disp('Fixed point iteration converged at every time step')
end

Total_time = toc;

disp(['Total time = ' num2str(Total_time)])

filename = ['data/particle_method_DG_2d_n',num2str(n),'_dv',num2str(dv),'_dt',...
    num2str(dt),'_t0_',num2str(t0),'_time_',num2str(time),'.mat'];
save(filename)
end
