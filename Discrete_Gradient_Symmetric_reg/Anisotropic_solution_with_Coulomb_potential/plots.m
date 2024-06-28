clc 
clear all
close all

data_1 = load('data/particle_2d_n40_dv0.5_dt_0.05_time_20.mat');
data_2 = load('data/particle_2d_n45_dv0.44444_dt_0.05_time_20.mat');
data_3 = load('data/particle_2d_n50_dv0.4_dt_0.05_time_20.mat');
data_4 = load('data/particle_2d_n55_dv0.36364_dt_0.05_time_20.mat');
data_5 = load('data/particle_2d_n60_dv0.33333_dt_0.05_time_20.mat');


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


%Weights, Initial Velocities, Initial energy ect...

vx_1 = data_1.vx;
vy_1 = data_1.vy;
Vx0_1 = vx_1(:);
Vy0_1 = vy_1(:);
W_1 = data_1.W;

E0_1 = sum(W_1.*(Vx0_1.^2+Vy0_1.^2));
Mx0_1 = sum(W_1.*Vx0_1);
My0_1 = sum(W_1.*Vy0_1);


vx_2 = data_2.vx;
vy_2 = data_2.vy;
Vx0_2 = vx_2(:);
Vy0_2 = vy_2(:);
W_2 = data_2.W;

E0_2 = sum(W_2.*(Vx0_2.^2+Vy0_2.^2));
Mx0_2 = sum(W_2.*Vx0_2);
My0_2 = sum(W_2.*Vy0_2);

vx_3 = data_3.vx;
vy_3 = data_3.vy;
Vx0_3 = vx_3(:);
Vy0_3 = vy_3(:);
W_3 = data_3.W;

E0_3 = sum(W_3.*(Vx0_3.^2+Vy0_3.^2));
Mx0_3 = sum(W_3.*Vx0_3);
My0_3 = sum(W_3.*Vy0_3);

vx_4 = data_4.vx;
vy_4 = data_4.vy;
Vx0_4 = vx_4(:);
Vy0_4 = vy_4(:);
W_4 = data_4.W;

E0_4 = sum(W_4.*(Vx0_4.^2+Vy0_4.^2));
Mx0_4 = sum(W_4.*Vx0_4);
My0_4 = sum(W_4.*Vy0_4);

vx_5 = data_5.vx;
vy_5 = data_5.vy;
Vx0_5 = vx_5(:);
Vy0_5 = vy_5(:);
W_5 = data_5.W;

E0_5 = sum(W_5.*(Vx0_5.^2+Vy0_5.^2));
Mx0_5 = sum(W_5.*Vx0_5);
My0_5 = sum(W_5.*Vy0_5);

figure 
plot(error_list_1(:,1),error_list_1(:,5)-E0_1,'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),......
    'DisplayName',['M = ',num2str(n_1)])
hold on 
plot(error_list_2(:,1),error_list_2(:,5)-E0_2,'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),......
    'DisplayName',['M = ',num2str(n_2)])
hold on 
plot(error_list_3(:,1),error_list_3(:,5)-E0_3,'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),......
    'DisplayName',['M = ',num2str(n_3)])
hold on 
plot(error_list_4(:,1),error_list_4(:,5)-E0_4,'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),......
    'DisplayName',['M = ',num2str(n_4)])
hold on 
plot(error_list_5(:,1),error_list_5(:,5)-E0_5,'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),......
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('Energy','Interpreter','latex')
title('Kinetic Energy Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_Energy.pdf','-dpdf','-r0')

%Entropy

figure 
plot(error_list_1(:,1),error_list_1(:,6),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,6),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,6),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,6),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,6),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^{\varepsilon}_{L}(f^N)$','Interpreter','latex')
title('Entropy Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_Entropy.pdf','-dpdf','-r0')

%Fisher Information

figure 
plot(error_list_1(:,1),error_list_1(:,8),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,8),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,8),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,8),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,8),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$F^{\varepsilon}(f^N)$','Interpreter','latex')
title('Fisher Information Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_Fisher.pdf','-dpdf','-r0')

%Dissipation term

figure 
plot(error_list_1(:,1),error_list_1(:,9),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,9),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,9),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,9),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,9),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$-D^{\varepsilon}(f^N)$','Interpreter','latex')
title('Dissipation Term Evolution','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_Dissipation.pdf','-dpdf','-r0')

figure
plot(error_list_1(:,1),error_list_1(:,3)-Mx0_1,'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,3)-Mx0_2,'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,3)-Mx0_3,'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,3)-Mx0_4,'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,3)-Mx0_5,'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('Momentum','Interpreter','latex')
title('Momentum Evolution $x$ Direction','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_momentum_x.pdf','-dpdf','-r0')

figure
plot(error_list_1(:,1),error_list_1(:,4)-My0_1,'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
plot(error_list_2(:,1),error_list_2(:,4)-My0_2,'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
plot(error_list_3(:,1),error_list_3(:,4)-My0_3,'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
plot(error_list_4(:,1),error_list_4(:,4)-My0_4,'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
plot(error_list_5(:,1),error_list_5(:,4)-My0_5,'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off
xlabel('$t$','Interpreter','latex')
ylabel('Momentum','Interpreter','latex')
title('Momentum Evolution $y$ Direction','Interpreter','latex')
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_momentum_y.pdf','-dpdf','-r0')

%Fixed Point Iteration Study

figure
semilogy(error_list_1(2:Nt_1,1),error_list_1(2:Nt_1,10),'-o','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_1,10)),...
    'DisplayName',['M = ',num2str(n_1)])
hold on
semilogy(error_list_2(2:Nt_2,1),error_list_2(2:Nt_2,10),'-*','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_2,10)),...
    'DisplayName',['M = ',num2str(n_2)])
hold on
semilogy(error_list_3(2:Nt_3,1),error_list_3(2:Nt_3,10),'-square','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_3,10)),...
    'DisplayName',['M = ',num2str(n_3)])
hold on
semilogy(error_list_4(2:Nt_4,1),error_list_4(2:Nt_4,10),'-+','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_4,10)),...
    'DisplayName',['M = ',num2str(n_4)])
hold on
semilogy(error_list_5(2:Nt_5,1),error_list_5(2:Nt_5,10),'-diamond','LineWidth',2,'MarkerIndices',...
    round(linspace(1,Nt_5,10)),...
    'DisplayName',['M = ',num2str(n_5)])
hold off 
hl = legend('show');
set(hl,'Location','northeast')
set(hl,'Interpreter','latex')
xlabel('t','Interpreter','latex')
ylabel('Iterations','Interpreter','latex')
title('Fixed Point Iterations','Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'plots/Coulomb_DG_eq_iter.pdf','-dpdf','-r0')

[max_1,max_1_ind] = max(error_list_1(:,10));
disp(['The maximum number of iterations for n = ',num2str(n_1),...
    ' is ',num2str(max_1),...
    ' and it occurs at time ',num2str(error_list_1(max_1_ind,1))])
avg_1 = mean(error_list_1(:,10));
disp(['The average number of iterations for n = ',num2str(n_1),...
    ' is ',num2str(avg_1)])


[max_2,max_2_ind] = max(error_list_2(:,10));
disp(['The maximum number of iterations for n = ',num2str(n_2),...
    ' is ',num2str(max_2),...
    ' and it occurs at time ',num2str(error_list_2(max_2_ind,1))])
avg_2 = mean(error_list_2(:,10));
disp(['The average number of iterations for n = ',num2str(n_2),...
    ' is ',num2str(avg_2)])

[max_3,max_3_ind] = max(error_list_3(:,10));
disp(['The maximum number of iterations for n = ',num2str(n_3),...
    ' is ',num2str(max_3),...
    ' and it occurs at time ',num2str(error_list_3(max_3_ind,1))])
avg_3 = mean(error_list_3(:,10));
disp(['The average number of iterations for n = ',num2str(n_3),...
    ' is ',num2str(avg_3)])

[max_4,max_4_ind] = max(error_list_4(:,10));
disp(['The maximum number of iterations for n = ',num2str(n_4),...
    ' is ',num2str(max_4),...
    ' and it occurs at time ',num2str(error_list_4(max_4_ind,1))])
avg_4 = mean(error_list_4(:,10));
disp(['The average number of iterations for n = ',num2str(n_4),...
    ' is ',num2str(avg_4)])

[max_5,max_5_ind] = max(error_list_5(:,10));
disp(['The maximum number of iterations for n = ',num2str(n_5),...
    ' is ',num2str(max_5),...
    ' and it occurs at time ',num2str(error_list_5(max_5_ind,1))])
avg_5 = mean(error_list_5(:,10));
disp(['The average number of iterations for n = ',num2str(n_5),...
    ' is ',num2str(avg_5)])




