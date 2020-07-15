% The following code computes the transient response of the TLSC SAPF 
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

%% Closed-loop performance indexes

% PI - Pole Placement

load('PIPP_TimeResponse_20200102.mat')

tPIPP = PIPPTimeResponse20200102(:,1);
idPIPP = PIPPTimeResponse20200102(:,2);
idPIPPRef = PIPPTimeResponse20200102(:,3);

iqPIPP = PIPPTimeResponse20200102(:,4);
iqPIPPRef = PIPPTimeResponse20200102(:,5);

i0PIPP = PIPPTimeResponse20200102(:,6);
i0PIPPRef = PIPPTimeResponse20200102(:,7);

vdcPIPP = PIPPTimeResponse20200102(:,8);
vdcPIPPRef = PIPPTimeResponse20200102(:,9);


load('PIOpt_TimeResponse_20200102.mat')

tPIOpt = PIOptTimeResponse20200102(:,1);
idPIOpt = PIOptTimeResponse20200102(:,2);
idPIOptRef = PIOptTimeResponse20200102(:,3);

iqPIOpt = PIOptTimeResponse20200102(:,4);
iqPIOptRef = PIOptTimeResponse20200102(:,5);

i0PIOpt = PIOptTimeResponse20200102(:,6);
i0PIOptRef = PIOptTimeResponse20200102(:,7);

vdcPIOpt = PIOptTimeResponse20200102(:,8);
vdcPIOptRef = PIOptTimeResponse20200102(:,9);

load('PIRobust_TimeResponse_20200102.mat')

tPIRobust = PIRobustTimeResponse20200102(:,1);
idPIRobust = PIRobustTimeResponse20200102(:,2);
idPIRobustRef = PIRobustTimeResponse20200102(:,3);

iqPIRobust = PIRobustTimeResponse20200102(:,4);
iqPIRobustRef = PIRobustTimeResponse20200102(:,5);

i0PIRobust = PIRobustTimeResponse20200102(:,6);
i0PIRobustRef = PIRobustTimeResponse20200102(:,7);

vdcPIRobust = PIRobustTimeResponse20200102(:,8);
vdcPIRobustRef = PIRobustTimeResponse20200102(:,9);

figure(1)
plot(tPIPP,idPIPP,'-o','LineWidth', 1); hold all;
plot(tPIOpt,idPIOpt,':d','LineWidth', 1);
plot(tPIRobust,idPIRobust,'--','LineWidth', 4);
plot(tPIPP,idPIPPRef,'LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{s}^{d}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-15 20])
set(gca,'FontSize',32,'FontName','Times')
legend_id = legend('PP-PI','$H\infty$-PI','$H\infty$-uPI','$i_{s_{REF}}^{d}$','show');
set(legend_id,'Interpreter','latex','Location','best');
hold off

figure(2)
plot(tPIPP,iqPIPP,'-o','LineWidth', 1); hold all;
plot(tPIOpt,iqPIOpt,':d','LineWidth', 1);
plot(tPIRobust,iqPIRobust,'--','LineWidth', 4);
plot(tPIPP,iqPIPPRef,'LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{s}^{q}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-6 5])
set(gca,'FontSize',32,'FontName','Times')
legend_iq = legend('PP-PI','$H\infty$-PI','$H\infty$-uPI','$i_{s_{REF}}^{q}$','show');
set(legend_iq,'Interpreter','latex','Location','best');
hold off

figure(3)
plot(tPIPP,i0PIPP,'-o','LineWidth', 1); hold all;
plot(tPIOpt,i0PIOpt,':d','LineWidth', 1);
plot(tPIRobust,i0PIRobust,'--','LineWidth', 4);
plot(tPIPP,i0PIPPRef,'LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{s}^{0}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-2 2.4])
set(gca,'FontSize',32,'FontName','Times')
legend_i0 = legend('PP-PI','$H\infty$-PI','$H\infty$-uPI','$i_{s_{REF}}^{0}$','show');
set(legend_i0,'Interpreter','latex','Location','best');
hold off

figure(4)
plot(tPIPP,vdcPIPP,'-o','LineWidth', 1); hold all;
plot(tPIOpt,vdcPIOpt,':d','LineWidth', 1);
plot(tPIRobust,vdcPIRobust,'--','LineWidth', 4);
plot(tPIPP,vdcPIPPRef,'LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$v_{dc}$ (V)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([580 625])
set(gca,'FontSize',32,'FontName','Times')
legend_vdc = legend('PP-PI','$H\infty$-PI','$H\infty$-uPI','$v_{dc_{REF}}$','show');
set(legend_vdc,'Interpreter','latex','Location','best');
hold off