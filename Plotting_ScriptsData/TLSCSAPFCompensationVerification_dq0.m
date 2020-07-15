% The following code plots figures related to the dq0 response 
% reported in "Robust Control for Shunt Active Power Filters: 
% A Dynamical Model-based Approach with VerifiedControllability", 
% by Jorge-Humberto Urrea-Quintero, Jesús M. López-Lezama
% Nicolás Muñoz-Galeano, Energies, 2020.

% It is a .m file, and therefore it is meant to be used on Matlab.
% To plot the comparisons, it is required to import the data on the .txt
% files as .mat files.

% Code written by: Jorge-Humberto Urrea-Quintero.
% 2020

%%
clear all;
close all;
clc;
%% Verification of the SAPF compensation dq0

%% id components

load('TLSCCurrentsdq0RobustPI_20200118.mat')

tidq0 = TLSCCurrentsdq0RobustPI20200118(:,1);
idRed = TLSCCurrentsdq0RobustPI20200118(:,3);
iqRed = TLSCCurrentsdq0RobustPI20200118(:,6);
i0Red = TLSCCurrentsdq0RobustPI20200118(:,10);
inRed = TLSCCurrentsdq0RobustPI20200118(:,13);
iLd = TLSCCurrentsdq0RobustPI20200118(:,2);
iLq = TLSCCurrentsdq0RobustPI20200118(:,5);
iL0 = TLSCCurrentsdq0RobustPI20200118(:,8);
iLn = TLSCCurrentsdq0RobustPI20200118(:,11);
ipd = TLSCCurrentsdq0RobustPI20200118(:,4);
ipq = TLSCCurrentsdq0RobustPI20200118(:,7);
ip0 = TLSCCurrentsdq0RobustPI20200118(:,9);
ipn = TLSCCurrentsdq0RobustPI20200118(:,12);

% id
figure(1)
plot(tidq0,idRed,'-','LineWidth', 1); hold all;
plot(tidq0,iLd,':','LineWidth', 1);
plot(tidq0,ipd,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$d$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-10 20])
set(gca,'FontSize',32,'FontName','Times')
legend_iRedabc = legend('$i_{grid}^{d}$','$i_{Load}^{d}$','$i_{S}^{d}$');
set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

%iq
figure(2)
plot(tidq0,iqRed,'-','LineWidth', 1); hold all;
plot(tidq0,iLq,':','LineWidth', 1);
plot(tidq0,ipq,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$q$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-5 5])
set(gca,'FontSize',32,'FontName','Times')
legend_iRedabc = legend('$i_{grid}^{q}$','$i_{Load}^{q}$','$i_{S}^{q}$');
set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(5)
plot(tidq0,iqRed,'-','LineWidth', 1); hold all;
plot(tidq0,iLq,':','LineWidth', 1);
plot(tidq0,ipq,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$q$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-2 2])
set(gca,'FontSize',32,'FontName','Times')
% legend_iRedabc = legend('$i_{grid}^{q}$','$i_{Load}^{q}$','$i_{S}^{q}$');
% set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(6)
plot(tidq0,iqRed,'-','LineWidth', 1); hold all;
plot(tidq0,iLq,':','LineWidth', 1);
plot(tidq0,ipq,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$q$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-1.5 1.5])
set(gca,'FontSize',32,'FontName','Times')
% legend_iRedabc = legend('$i_{grid}^{q}$','$i_{Load}^{q}$','$i_{S}^{q}$');
% set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

%i0
figure(3)
plot(tidq0,i0Red,'-','LineWidth', 1); hold all;
plot(tidq0,iL0,':','LineWidth', 1);
plot(tidq0,ip0,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-2 2])
set(gca,'FontSize',32,'FontName','Times')
legend_iRedabc = legend('$i_{grid}^{0}$','$i_{Load}^{0}$','$i_{S}^{0}$');
set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% Zoom
figure(7)
plot(tidq0,i0Red,'-','LineWidth', 1); hold all;
plot(tidq0,iL0,':','LineWidth', 1);
plot(tidq0,ip0,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0$ component (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([1.4 1.5])
ylim([-2 2])
set(gca,'FontSize',32,'FontName','Times')
% legend_iRedabc = legend('$i_{grid}^{0}$','$i_{Load}^{0}$','$i_{S}^{0}$');
% set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

%in
figure(4)
plot(tidq0,inRed,'-','LineWidth', 1); hold all;
plot(tidq0,iLn,':','LineWidth', 1);
plot(tidq0,ipn,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('neutral current (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-6 6])
set(gca,'FontSize',32,'FontName','Times')
legend_iRedabc = legend('$i_{grid}^{n}$','$i_{Load}^{n}$','$i_{S}^{n}$');
set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

