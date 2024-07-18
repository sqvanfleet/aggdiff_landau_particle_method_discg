clear all
close all
clc


data_60 = load('Data_m_2/particle_2d_n60_dx0.26667_epsilon0.04673_dt0.01_tmax3_m_2.mat');
data_70 = load('Data_m_2/particle_2d_n70_dx0.22857_epsilon0.034438_dt0.01_tmax3_m_2.mat');
data_80 = load('Data_m_2/particle_2d_n80_dx0.2_epsilon0.026437_dt0.01_tmax3_m_2.mat');
data_90 = load('Data_m_2/particle_2d_n90_dx0.17778_epsilon0.020938_dt0.01_tmax3_m_2.mat');
data_100 = load('Data_m_2/particle_2d_n100_dx0.16_epsilon0.016996_dt0.01_tmax3_m_2.mat');

n_60 = data_60.n;
n_70 = data_70.n;
n_80 = data_80.n;
n_90 = data_90.n;
n_100 = data_100.n;

Nt_60 = data_60.Nt;
Nt_70 = data_70.Nt;
Nt_80 = data_80.Nt;
Nt_90 = data_90.Nt;
Nt_100 = data_100.Nt;

error_list_60 = data_60.error_list;
error_list_70 = data_70.error_list;
error_list_80 = data_80.error_list;
error_list_90 = data_90.error_list;
error_list_100 = data_100.error_list;

%Error Evolution L2


figure 
plot(error_list_60(:,1),error_list_60(:,4),'-o','MarkerIndices',...
    round(linspace(1,Nt_60,10)),...
    'DisplayName',['n = ',num2str(n_60)])
hold on 
plot(error_list_70(:,1),error_list_70(:,4),'-*','MarkerIndices',...
    round(linspace(1,Nt_70,10)),...
    'DisplayName',['n = ',num2str(n_70)])
hold on 
plot(error_list_80(:,1),error_list_80(:,4),'-square','MarkerIndices',...
    round(linspace(1,Nt_80,10)),...
    'DisplayName',['n = ',num2str(n_80)])
hold on 
plot(error_list_90(:,1),error_list_90(:,4),'-+','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold on 
plot(error_list_100(:,1),error_list_100(:,4),'-diamond','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold off
legend 
xlabel('t')
ylabel('f')
title('L2 Error Evolution')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_porous_medium_L2_error.pdf','-dpdf','-r0')

%Error Evolution L1


figure 
plot(error_list_60(:,1),error_list_60(:,3),'-o','MarkerIndices',...
    round(linspace(1,Nt_60,10)),...
    'DisplayName',['n = ',num2str(n_60)])
hold on 
plot(error_list_70(:,1),error_list_70(:,3),'-*','MarkerIndices',...
    round(linspace(1,Nt_70,10)),...
    'DisplayName',['n = ',num2str(n_70)])
hold on 
plot(error_list_80(:,1),error_list_80(:,3),'-square','MarkerIndices',...
    round(linspace(1,Nt_80,10)),...
    'DisplayName',['n = ',num2str(n_80)])
hold on 
plot(error_list_90(:,1),error_list_90(:,3),'-+','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold on 
plot(error_list_100(:,1),error_list_100(:,3),'-diamond','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold off
legend 
xlabel('t')
ylabel('f')
title('L1 Error Evolution')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_porous_medium_L1_error.pdf','-dpdf','-r0')

%Error Evolution Linf

figure 
plot(error_list_60(:,1),error_list_60(:,2),'-o','MarkerIndices',...
    round(linspace(1,Nt_60,10)),...
    'DisplayName',['n = ',num2str(n_60)])
hold on 
plot(error_list_70(:,1),error_list_70(:,2),'-*','MarkerIndices',...
    round(linspace(1,Nt_70,10)),...
    'DisplayName',['n = ',num2str(n_70)])
hold on 
plot(error_list_80(:,1),error_list_80(:,2),'-square','MarkerIndices',...
    round(linspace(1,Nt_80,10)),...
    'DisplayName',['n = ',num2str(n_80)])
hold on 
plot(error_list_90(:,1),error_list_90(:,2),'-+','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold on 
plot(error_list_100(:,1),error_list_100(:,2),'-diamond','MarkerIndices',...
    round(linspace(1,Nt_90,10)),...
    'DisplayName',['n = ',num2str(n_90)])
hold off
legend 
xlabel('t')
ylabel('f')
title('Linf Error Evolution')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_porous_medium_Linf_error.pdf','-dpdf','-r0')



%Mass evolution 

% figure
% plot(error_list_60(:,1),error_list_60(:,5),'DisplayName','Mass')
% legend
% xlabel('t')
% ylabel('$\rho$','Interpreter','latex')
% title('Mass')


% %L2 error evolution comparison
% 
% figure 
% plot(error_list_30(:,1),error_list_30(:,4),'-o','MarkerIndices',...
%     round(linspace(1,Nt(1),10)),...
%     'DisplayName',['n = ',num2str(n(1))])
% hold on
% plot(error_list_60(:,1),error_list_60(:,4),'-*','MarkerIndices',...
%     round(linspace(1,Nt(2),10)),...
%     'DisplayName',['n = ',num2str(n(2))])
% hold on
% plot(error_list_120(:,1),error_list_120(:,4),'-square','MarkerIndices',...
%     round(linspace(1,Nt(3),10)),...
%     'DisplayName',['n = ',num2str(n(3))])
% hold on
% plot(error_list_240(:,1),error_list_240(:,4),'-^','MarkerIndices',...
%     round(linspace(1,Nt(2),10)),...
%     'DisplayName',['n = ',num2str(n(4))])
% hold off 
% legend
% xlabel('t','Interpreter','latex')
% ylabel('$L_2 Error$','Interpreter','latex')
% title('Error Evolution')
% f = gcf;
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f,'1D_wave_order_of_accuracy.pdf','-dpdf','-r0')

%Entropy Evolution

figure 
semilogy(error_list_60(:,1),error_list_60(:,6),...
    '-o','MarkerIndices',round(linspace(1,Nt_60,10)),...
    'DisplayName',['Entropy n = ',num2str(n_60)])
hold on
semilogy(error_list_70(:,1),error_list_70(:,6),...
    '-*','MarkerIndices',round(linspace(1,Nt_70,10)),...
    'DisplayName',['Entropy n = ',num2str(n_70)])
hold on
semilogy(error_list_80(:,1),error_list_80(:,6),...
    '-square','MarkerIndices',round(linspace(1,Nt_80,10)),...
    'DisplayName',['Entropy n = ',num2str(n_80)])
hold on
semilogy(error_list_90(:,1),error_list_90(:,6),...
    '-+','MarkerIndices',round(linspace(1,Nt_90,10)),...
    'DisplayName',['Entropy n = ',num2str(n_90)])
hold on
semilogy(error_list_100(:,1),error_list_100(:,6),...
    '-diamond','MarkerIndices',round(linspace(1,Nt_100,10)),...
    'DisplayName',['Entropy n = ',num2str(n_100)])
hold off
legend
xlabel('t')
ylabel('E(f)')
title('Entropy')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_porous_medium_Entropy.pdf','-dpdf','-r0')
% f = gcf;
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f,'1D_wave_entropy.pdf','-dpdf','-r0')




 
%error_ratio = zeros(3,1);

%L2 error list at t = 3
L2_error = [error_list_60(Nt_60,4),...
    error_list_70(Nt_70,4),error_list_80(Nt_80,4),...
    error_list_90(Nt_90,4),error_list_100(Nt_100,4)];

L1_error = [error_list_60(Nt_60,3),...
    error_list_70(Nt_70,3),error_list_80(Nt_80,3),...
    error_list_90(Nt_90,3),error_list_100(Nt_100,3)];

Linf_error = [error_list_60(Nt_60,2),...
    error_list_70(Nt_70,2),error_list_80(Nt_80,2),...
    error_list_90(Nt_90,2),error_list_100(Nt_100,2)];

h_60 = data_60.dx; 
h_70 = data_70.dx;
h_80 = data_80.dx;
h_90 = data_90.dx;
h_100 = data_100.dx;

h_list = [h_60,h_70,h_80,h_90,h_100];

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
xlabel('$\log{dx}$','Interpreter','latex')
ylabel('$\log{\mbox{error}}$','Interpreter','latex')
title('Rate of Convergence $L^2$','Interpreter','latex')
hl = legend('show','Location','northwest');
set(hl,'Interpreter','latex')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'1D_porous_medium_ROC.pdf','-dpdf','-r0')


%Order of convergence at time t = 3
% for i = 1:length(error_ratio)
%     error_ratio(i) = (1/log(2))*log(L2_error(i)/L2_error(i+1));
% end
% disp(error_ratio);
% disp(L2_error');








