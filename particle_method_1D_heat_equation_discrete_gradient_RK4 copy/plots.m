clear all
close all
clc

data_1 = load('particle_2d_n60_dx0.5_epsilon0.16223_dt0.01_tmax3.mat');
data_2 = load('particle_2d_n70_dx0.42857_epsilon0.11956_dt0.01_tmax3.mat');
data_3 = load('particle_2d_n80_dx0.375_epsilon0.091783_dt0.01_tmax3.mat');
data_4 = load('particle_2d_n90_dx0.33333_epsilon0.072691_dt0.01_tmax3.mat');
data_5 = load('particle_2d_n100_dx0.3_epsilon0.059004_dt0.01_tmax3.mat');


error_list_1 = data_1.error_list;
error_list_2 = data_2.error_list;
error_list_3 = data_3.error_list;
error_list_4 = data_4.error_list;
error_list_5 = data_5.error_list;

n_1 = data_1.n; n_2 = data_2.n; n_3 = data_3.n;
n_4 = data_4.n; n_5 = data_5.n;

dx_1 = data_1.dx; dx_2 = data_2.dx; dx_3 = data_3.dx;
dx_4 = data_4.dx; dx_5 = data_5.dx;

Nt_1 = data_1.Nt; Nt_2 = data_2.Nt; Nt_3 = data_3.Nt;
Nt_4 = data_4.Nt; Nt_5 = data_5.Nt;

%L2 error evolution comparison

figure 
plot(error_list_1(:,1),error_list_1(:,4),'-o','MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['n = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,4),'-*','MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['n = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,4),'-square','MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['n = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,4),'-+','MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['n = ',num2str(n_4)])
plot(error_list_5(:,1),error_list_5(:,4),'-diamond','MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['n = ',num2str(n_5)])
hold off 
legend
xlabel('$t$','Interpreter','latex')
ylabel('$L^2$ Error','Interpreter','latex')
title('$L^2$ Error Evolution','Interpreter','latex')
f = gcf;
fontsize(legend,16,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_heat_eq_rk4_error.pdf','-dpdf','-r0')

%Entropy Evolution

figure 
plot(error_list_1(:,1),error_list_1(:,6),...
    '-o','MarkerIndices',round(linspace(1,Nt_1,10)),...
    'DisplayName',['Entropy n = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,6),...
    '-*','MarkerIndices',round(linspace(1,Nt_2,10)),...
    'DisplayName',['Entropy n = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,6),...
    '-square','MarkerIndices',round(linspace(1,Nt_3,10)),...
    'DisplayName',['Entropy n = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,6),...
    '-+','MarkerIndices',round(linspace(1,Nt_4,10)),...
    'DisplayName',['Entropy n = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,6),...
    '-diamond','MarkerIndices',round(linspace(1,Nt_5,10)),...
    'DisplayName',['Entropy n = ',num2str(n_5)])
hold off
legend
xlabel('$t$','Interpreter','latex')
ylabel('$F^{\varepsilon}_{A}(f^N)$','Interpreter','latex')
title('\textbf{Energy Evolution}','Interpreter','latex')
f = gcf;
fontsize(legend,14,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_heat_eq_rk4_entropy.pdf','-dpdf','-r0')

%L2 error list at t = 3
L2_error = [error_list_1(Nt_1,4),...
    error_list_2(Nt_2,4),error_list_3(Nt_3,4),...
    error_list_4(Nt_4,4),error_list_5(Nt_5,4)];

L1_error = [error_list_1(Nt_1,3),...
    error_list_2(Nt_2,3),error_list_3(Nt_3,3),...
    error_list_4(Nt_4,3),error_list_5(Nt_5,3)];

Linf_error = [error_list_1(Nt_1,2),...
    error_list_2(Nt_2,2),error_list_3(Nt_3,2),...
    error_list_4(Nt_4,2),error_list_5(Nt_5,2)];

h_1 = data_1.dx; 
h_2 = data_2.dx;
h_3 = data_3.dx;
h_4 = data_4.dx;
h_5 = data_5.dx;


h_list = [h_1,h_2,h_3,h_4,h_5];

x = ones(length(h_list),2);
x(:,2) = log(h_list);

L2_roc = x\log(L2_error'); %slope using least squares
L1_roc = x\log(L1_error'); 
Linf_roc = x\log(Linf_error');

figure 
plot(-log(h_list),-L2_roc(2)*log(h_list) - L2_roc(1),...
    '-o','color','r','DisplayName',['$L^2$ error slope $=$ ',num2str(L2_roc(2))])
hold on 
plot(-log(h_list),-L1_roc(2)*log(h_list) - L1_roc(1),...
    '-square','color','b','DisplayName',['$L^1$ error slope $=$ ',num2str(L1_roc(2))])
hold on 
plot(-log(h_list),-Linf_roc(2)*log(h_list) - Linf_roc(1),...
    '-diamond','color','k','DisplayName',['$L^{\infty}$ error slope $=$ ',num2str(Linf_roc(2))])
hold on 
plot(-log(h_list),-2*log(h_list)+3,"-.",'DisplayName','Reference Slope')
xlabel('$\log{(h)}$','Interpreter','latex')
ylabel('$\log{(\mbox{error})}$','Interpreter','latex')
title('\textbf{Rate of Convergence}','Interpreter','latex')
hl = legend('show','Location','northwest');
fontsize(hl,14,"points")
set(hl,'Interpreter','latex')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_heat_eq_rk4_ROC.pdf','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure
semilogy(error_list_1(:,1),error_list_1(:,7),'-o','MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['n = ',num2str(n_1)])
hold on
semilogy(error_list_2(:,1),error_list_2(:,7),'-*','MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['n = ',num2str(n_2)])
hold on
semilogy(error_list_3(:,1),error_list_3(:,7),'-square','MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['n = ',num2str(n_3)])
hold on
semilogy(error_list_4(:,1),error_list_4(:,7),'-+','MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['n = ',num2str(n_4)])
hold on
semilogy(error_list_5(:,1),error_list_5(:,7),'-diamond','MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['n = ',num2str(n_5)])
hold off 
legend
xlabel('$t$','Interpreter','latex')
ylabel('$\mbox{Iterations}$','Interpreter','latex')
title('\textbf{Fixed Point Iterations}','Interpreter','latex')
f = gcf;
fontsize(legend,14,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_heat_eq_rk4_iter.pdf','-dpdf','-r0')

[max_1,max_1_ind] = max(error_list_1(:,7));
disp(['The maximum number of iterations for n = ',num2str(n_1),...
    ' is ',num2str(max_1),...
    ' and it occurs at time ',num2str(error_list_1(max_1_ind,1))])
avg_1 = mean(error_list_1(:,7));
disp(['The average number of iterations for n = ',num2str(n_1),...
    ' is ',num2str(avg_1)])


[max_2,max_2_ind] = max(error_list_2(:,7));
disp(['The maximum number of iterations for n = ',num2str(n_2),...
    ' is ',num2str(max_2),...
    ' and it occurs at time ',num2str(error_list_2(max_2_ind,1))])
avg_2 = mean(error_list_2(:,7));
disp(['The average number of iterations for n = ',num2str(n_2),...
    ' is ',num2str(avg_2)])

[max_3,max_3_ind] = max(error_list_3(:,7));
disp(['The maximum number of iterations for n = ',num2str(n_3),...
    ' is ',num2str(max_3),...
    ' and it occurs at time ',num2str(error_list_3(max_3_ind,1))])
avg_3 = mean(error_list_3(:,7));
disp(['The average number of iterations for n = ',num2str(n_3),...
    ' is ',num2str(avg_3)])

[max_4,max_4_ind] = max(error_list_4(:,7));
disp(['The maximum number of iterations for n = ',num2str(n_4),...
    ' is ',num2str(max_4),...
    ' and it occurs at time ',num2str(error_list_4(max_4_ind,1))])
avg_4 = mean(error_list_4(:,7));
disp(['The average number of iterations for n = ',num2str(n_4),...
    ' is ',num2str(avg_4)])

[max_5,max_5_ind] = max(error_list_5(:,7));
disp(['The maximum number of iterations for n = ',num2str(n_5),...
    ' is ',num2str(max_5),...
    ' and it occurs at time ',num2str(error_list_5(max_5_ind,1))])
avg_5 = mean(error_list_5(:,7));
disp(['The average number of iterations for n = ',num2str(n_5),...
    ' is ',num2str(avg_5)])



% figure
% scatter(error_list_30(:,1),error_list_30(:,7),...
%     'DisplayName',['n = ',num2str(n(1))])
% hold on
% scatter(error_list_60(:,1),error_list_60(:,7),...
%     'DisplayName',['n = ',num2str(n(2))])
% hold on 
% scatter(error_list_120(:,1),error_list_120(:,7),...
%     'DisplayName',['n = ',num2str(n(1))])
% hold on
% scatter(error_list_240(:,1),error_list_240(:,7),...
%     'DisplayName',['n = ',num2str(n(1))])
% hold off 
% legend
% xlabel('t','Interpreter','latex')
% ylabel('$\mbox{Iterations}$','Interpreter','latex')
% title('Number of iterations')
% f = gcf;
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f,'1D_heat_eq_iter.pdf','-dpdf','-r0')









