clc 
clear all
close all

data_1 = load('data/particle_method_DG_2d_n40_dv0.2_dt0.00125_t0_0_time_5.mat');
data_2 = load('data/particle_method_DG_2d_n45_dv0.17778_dt0.00125_t0_0_time_5.mat');
data_3 = load('data/particle_method_DG_2d_n50_dv0.16_dt0.00125_t0_0_time_5.mat');
data_4 = load('data/particle_method_DG_2d_n55_dv0.14545_dt0.00125_t0_0_time_5.mat');
data_5 = load('data/particle_method_DG_2d_n60_dv0.13333_dt0.00125_t0_0_time_5.mat');


Nt_1 = data_1.Nt;
Nt_2 = data_2.Nt;
Nt_3 = data_3.Nt;
Nt_4 = data_4.Nt;
Nt_5 = data_5.Nt;

n_1 = data_1.n;
n_2 = data_2.n;
n_3 = data_3.n;
n_4 = data_4.n;
n_5 = data_5.n;

error_list_1 = data_1.error_list;
error_list_2 = data_2.error_list;
error_list_3 = data_3.error_list;
error_list_4 = data_4.error_list;
error_list_5 = data_5.error_list;

figure 
plot(error_list_1(:,1),error_list_1(:,4),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['$M = $ ',num2str(n_1)])
hold on 
plot(error_list_2(:,1),error_list_2(:,4),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['$M = $ ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,4),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),......
    'DisplayName',['$M = $ ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,4),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),......
    'DisplayName',['$M = $ ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,4),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),......
    'DisplayName',['$M = $ ',num2str(n_5)])
xlabel('Time')
ylabel('L2 Error')
legend
title('L2 error evolution')




%Weights, Initial Velocities, Initial energy ect...

vx_1 = data_1.vx;
vy_1 = data_1.vy;
Vx0_1 = vx_1(:);
Vy0_1 = vy_1(:);
W_1 = data_1.W;

E0_1 = sum(W_1.*(Vx0_1.^2+Vy0_1.^2));


vx_2 = data_2.vx;
vy_2 = data_2.vy;
Vx0_2 = vx_2(:);
Vy0_2 = vy_2(:);
W_2 = data_2.W;

E0_2 = sum(W_2.*(Vx0_2.^2+Vy0_2.^2));

vx_3 = data_3.vx;
vy_3 = data_3.vy;
Vx0_3 = vx_3(:);
Vy0_3 = vy_3(:);
W_3 = data_3.W;

E0_3 = sum(W_3.*(Vx0_3.^2+Vy0_3.^2));

vx_4 = data_4.vx;
vy_4 = data_4.vy;
Vx0_4 = vx_4(:);
Vy0_4 = vy_4(:);
W_4 = data_4.W;

E0_4 = sum(W_4.*(Vx0_4.^2+Vy0_4.^2));

vx_5 = data_5.vx;
vy_5 = data_5.vy;
Vx0_5 = vx_5(:);
Vy0_5 = vy_5(:);
W_5 = data_5.W;

E0_5 = sum(W_5.*(Vx0_5.^2+Vy0_5.^2));

figure 
semilogy(error_list_1(:,1),abs(error_list_1(:,8)-E0_1),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),......
    'DisplayName',['$M = $ ',num2str(n_1)])
hold on 
semilogy(error_list_2(:,1),abs(error_list_2(:,8)-E0_2),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),......
    'DisplayName',['$M = $ ',num2str(n_2)])
hold on 
semilogy(error_list_3(:,1),abs(error_list_3(:,8)-E0_3),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),......
    'DisplayName',['$M = $ ',num2str(n_3)])
hold on 
semilogy(error_list_4(:,1),abs(error_list_4(:,8)-E0_4),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),......
    'DisplayName',['$M = $ ',num2str(n_4)])
hold on 
semilogy(error_list_5(:,1),abs(error_list_5(:,8)-E0_5),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),......
    'DisplayName',['$M = $ ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$|K-K^0|$','Interpreter','latex')
axis([0 5 10e-17 10e-15])
title('Kinetic Energy Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
ax = gca;
Ax.YLabel.Position(1) = -0.05;
f = gcf;
fontsize(f,28,"points")
fontsize(legend,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_Energy.pdf','-dpdf','-r0')

figure 
plot(error_list_1(:,1),error_list_1(:,9),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['$M = $ ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,9),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['$M = $ ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,9),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['$M = $ ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,9),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['$M = $ ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,9),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['$M = $ ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^{\varepsilon}_L$','Interpreter','latex')
axis tight
title('Energy Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,28,"points")
fontsize(legend,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_Entropy.pdf','-dpdf','-r0')

figure
plot(error_list_1(:,1),error_list_1(:,6),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['$M =$ ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,6),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['$M =$ ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,6),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['$M =$ ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,6),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['$M =$ ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,6),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['$M =$ ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('Momentum','Interpreter','latex')
title('Momentum Evolution $x$ Direction','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northwest')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_momentum_x.pdf','-dpdf','-r0')

figure
plot(error_list_1(:,1),error_list_1(:,7),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['$M =$ ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,7),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['$M =$ ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,7),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['$M =$ ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,7),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['$M =$ ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,7),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['$M =$ ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('Momentum','Interpreter','latex')
title('Momentum Evolution $y$ Direction','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northwest')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_momentum_y.pdf','-dpdf','-r0')




%%%% Log Log plot of the L1, L2, and Linfty error to display order of
%%%% convergence
P = 5; % Number of data sets
Linf_error = zeros(P,1);
L1_error = zeros(P,1);
L2_error = zeros(P,1);
h_list = zeros(P,1);

% The first column is h

h_list(1) = data_1.dv;
h_list(2) = data_2.dv;
h_list(3) = data_3.dv;
h_list(4) = data_4.dv;
h_list(5) = data_5.dv;


%the second column of error_list is the Linf error

Linf_error(1) = error_list_1(Nt_1,2);
Linf_error(2) = error_list_2(Nt_2,2);
Linf_error(3) = error_list_3(Nt_3,2);
Linf_error(4) = error_list_4(Nt_4,2);
Linf_error(5) = error_list_5(Nt_5,2);


%the third column of error_list is the L1 error
L1_error(1) = error_list_1(Nt_1,3);
L1_error(2) = error_list_2(Nt_2,3);
L1_error(3) = error_list_3(Nt_3,3);
L1_error(4) = error_list_4(Nt_4,3);
L1_error(5) = error_list_5(Nt_5,3);


%the fourth column of the error_list is the L2 error 
L2_error(1) = error_list_1(Nt_1,4);
L2_error(2) = error_list_2(Nt_2,4);
L2_error(3) = error_list_3(Nt_3,4);
L2_error(4) = error_list_4(Nt_4,4);
L2_error(5) = error_list_5(Nt_5,4);

x = ones(length(h_list),2);
x(:,2) = log(h_list);

Linf_roc = x\log(Linf_error); %slope using least squares fitting
L1_roc = x\log(L1_error); %slope using least squares fitting
L2_roc = x\log(L2_error); %slope using a least squares fitting


%Rate of convergence plot 
figure
plot(-log(h_list),-L2_roc(2)*log(h_list)-L2_roc(1),'-o','LineWidth',2,'DisplayName',...
    ['$L^2$ order $=$ ',num2str(L2_roc(2),4)])
hold on
plot(-log(h_list),-L1_roc(2)*log(h_list)-L1_roc(1),'-*','LineWidth',2,'DisplayName',...
    ['$L^1$ order $=$ ',num2str(L1_roc(2),4)])
hold on
plot(-log(h_list),-Linf_roc(2)*log(h_list)-Linf_roc(1)+.5,'-square','LineWidth',2,'DisplayName',...
    ['$L^{\infty}$ order $=$ ',num2str(Linf_roc(2),4)])
hold on
plot(-log(h_list),-2*log(h_list)+.4,'--','DisplayName','Reference Slope = 2')
xlabel('$\log{(h)}$','Interpreter','latex')
ylabel('$\log{\|f_{\varepsilon}^N - f\|}$','Interpreter','latex')
title('Rate of Convergence','Interpreter','latex')
hl = legend('show');
set(hl,'Location','southeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,28,"points")
fontsize(legend,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_ROC.pdf','-dpdf','-r0')

%Fixed point iteration study

figure
semilogy(error_list_1(2:Nt_1,1),error_list_1(2:Nt_1,13),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['$M = $ ',num2str(n_1)])
hold on
semilogy(error_list_2(2:Nt_2,1),error_list_2(2:Nt_2,13),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['$M = $ ',num2str(n_2)])
hold on
semilogy(error_list_3(2:Nt_3,1),error_list_3(2:Nt_3,13),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['$M = $ ',num2str(n_3)])
hold on
semilogy(error_list_4(2:Nt_4,1),error_list_4(2:Nt_4,13),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['$M = $ ',num2str(n_4)])
hold on
semilogy(error_list_5(2:Nt_5,1),error_list_5(2:Nt_5,13),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['$M = $ ',num2str(n_5)])
hold off 
hl = legend('show','Location','northeast');
set(hl,'Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('Iterations','Interpreter','latex')
axis tight
title('Fixed Point Iterations','Interpreter','latex')
f = gcf;
fontsize(f,28,"points")
fontsize(legend,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/BKW_DG_eq_iter.pdf','-dpdf','-r0')

[max_1,max_1_ind] = max(error_list_1(:,13));
disp(['The maximum number of iterations for M = ',num2str(n_1),...
    ' is ',num2str(max_1),...
    ' and it occurs at time ',num2str(error_list_1(max_1_ind,1))])
avg_1 = mean(error_list_1(:,13));
disp(['The average number of iterations for M = ',num2str(n_1),...
    ' is ',num2str(avg_1,4)])


[max_2,max_2_ind] = max(error_list_2(:,13));
disp(['The maximum number of iterations for M = ',num2str(n_2),...
    ' is ',num2str(max_2),...
    ' and it occurs at time ',num2str(error_list_2(max_2_ind,1))])
avg_2 = mean(error_list_2(:,13));
disp(['The average number of iterations for M = ',num2str(n_2),...
    ' is ',num2str(avg_2,4)])

[max_3,max_3_ind] = max(error_list_3(:,13));
disp(['The maximum number of iterations for M = ',num2str(n_3),...
    ' is ',num2str(max_3),...
    ' and it occurs at time ',num2str(error_list_3(max_3_ind,1))])
avg_3 = mean(error_list_3(:,13));
disp(['The average number of iterations for M = ',num2str(n_3),...
    ' is ',num2str(avg_3,4)])

[max_4,max_4_ind] = max(error_list_4(:,13));
disp(['The maximum number of iterations for M = ',num2str(n_4),...
    ' is ',num2str(max_4),...
    ' and it occurs at time ',num2str(error_list_4(max_4_ind,1))])
avg_4 = mean(error_list_4(:,13));
disp(['The average number of iterations for M = ',num2str(n_4),...
    ' is ',num2str(avg_4,4)])

[max_5,max_5_ind] = max(error_list_5(:,13));
disp(['The maximum number of iterations for M = ',num2str(n_5),...
    ' is ',num2str(max_5),...
    ' and it occurs at time ',num2str(error_list_5(max_5_ind,1))])
avg_5 = mean(error_list_5(:,13));
disp(['The average number of iterations for M = ',num2str(n_5),...
    ' is ',num2str(avg_5,4)])




