% The following code plots figures related to the abc response 
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
%% Verification of the SAPF compensation

%% Currents waveform

load('TLSCCurrentsRobustPI_20200117.mat')

tiabc = TLSCCurrentsRobustPI20200117(:,1);
iaRed = TLSCCurrentsRobustPI20200117(:,2);
ibRed = TLSCCurrentsRobustPI20200117(:,3);
icRed = TLSCCurrentsRobustPI20200117(:,4);
iLa = TLSCCurrentsRobustPI20200117(:,5);
iLb = TLSCCurrentsRobustPI20200117(:,6);
iLc = TLSCCurrentsRobustPI20200117(:,7);
ipa = TLSCCurrentsRobustPI20200117(:,8);
ipb = TLSCCurrentsRobustPI20200117(:,9);
ipc = TLSCCurrentsRobustPI20200117(:,10);

% igridabc
figure(1)
plot(tiabc,iaRed,'-','LineWidth', 1); hold all;
plot(tiabc,ibRed,':','LineWidth', 1);
plot(tiabc,icRed,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{grid}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_iRedabc = legend('$i_{grid}^{a}$','$i_{grid}^{b}$','$i_{grid}^{c}$');
set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(31)
plot(tiabc,iaRed,'-','LineWidth', 1); hold all;
plot(tiabc,ibRed,':','LineWidth', 1);
plot(tiabc,icRed,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{grid}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.6 0.7])
ylim([-10 10])
set(gca,'FontSize',32,'FontName','Times')
% legend_iRedabc = legend('$i_{grid}^{a}$','$i_{grid}^{b}$','$i_{grid}^{c}$');
% set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(32)
plot(tiabc,iaRed,'-','LineWidth', 1); hold all;
plot(tiabc,ibRed,':','LineWidth', 1);
plot(tiabc,icRed,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{grid}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-10 10])
set(gca,'FontSize',32,'FontName','Times')
% legend_iRedabc = legend('$i_{grid}^{a}$','$i_{grid}^{b}$','$i_{grid}^{c}$');
% set(legend_iRedabc,'Interpreter','latex','Location','best');
hold off

% iLoadabc
figure(2)
plot(tiabc,iLa,'-','LineWidth', 1); hold all;
plot(tiabc,iLb,':','LineWidth', 1);
plot(tiabc,iLc,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{Load}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.4])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_iLoadabc = legend('$i_{Load}^{a}$','$i_{Load}^{b}$','$i_{Load}^{c}$');
set(legend_iLoadabc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(33)
plot(tiabc,iLa,'-','LineWidth', 2); hold all;
plot(tiabc,iLb,':','LineWidth', 2);
plot(tiabc,iLc,'--','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{Load}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-10 10])
set(gca,'FontSize',32,'FontName','Times')
% legend_iLoadabc = legend('$i_{Load}^{a}$','$i_{Load}^{b}$','$i_{Load}^{c}$');
% set(legend_iLoadabc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(34)
plot(tiabc,iLa,'-','LineWidth', 2); hold all;
plot(tiabc,iLb,':','LineWidth', 2);
plot(tiabc,iLc,'--','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{Load}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-7 7])
set(gca,'FontSize',32,'FontName','Times')
% legend_iLoadabc = legend('$i_{Load}^{a}$','$i_{Load}^{b}$','$i_{Load}^{c}$');
% set(legend_iLoadabc,'Interpreter','latex','Location','best');
hold off

% iSabc
figure(3)
plot(tiabc,ipa,'-','LineWidth', 1); hold all;
plot(tiabc,ipb,':','LineWidth', 1);
plot(tiabc,ipc,'--','LineWidth', 1);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$i_{S}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_iSabc = legend('$i_{S}^{a}$','$i_{S}^{b}$','$i_{S}^{c}$');
set(legend_iSabc,'Interpreter','latex','Location','best');
hold off

% Zoom 1 
figure(35)
plot(tiabc,ipa,'-','LineWidth', 1); hold all;
plot(tiabc,ipb,':','LineWidth', 1);
plot(tiabc,ipc,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{S}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-5 5])
set(gca,'FontSize',32,'FontName','Times')
% legend_iSabc = legend('$i_{S}^{a}$','$i_{S}^{b}$','$i_{S}^{c}$');
% set(legend_iSabc,'Interpreter','latex','Location','best');
hold off

% Zoom 2 
figure(36)
plot(tiabc,ipa,'-','LineWidth', 1); hold all;
plot(tiabc,ipb,':','LineWidth', 1);
plot(tiabc,ipc,'--','LineWidth', 1);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$i_{S}^{abc}$ (A)', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-2 2])
set(gca,'FontSize',32,'FontName','Times')
% legend_iSabc = legend('$i_{S}^{a}$','$i_{S}^{b}$','$i_{S}^{c}$');
% set(legend_iSabc,'Interpreter','latex','Location','best');
hold off
%% Voltage load currents comparison
clear all

load('TLSCViLPhasesRobustPI_20200117.mat')

tViL = TLSCViLPhasesRobustPI20200117(:,1);
Vpcca = TLSCViLPhasesRobustPI20200117(:,2);
iLa = TLSCViLPhasesRobustPI20200117(:,3);
Vpccb = TLSCViLPhasesRobustPI20200117(:,4);
iLb = TLSCViLPhasesRobustPI20200117(:,5);
Vpccc = TLSCViLPhasesRobustPI20200117(:,6);
iLc = TLSCViLPhasesRobustPI20200117(:,7);

% iLa
figure(4)
plot(tViL,Vpcca,'-','LineWidth', 1); hold all;
plot(tViL,iLa,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Load}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{a}$','$i_{Load}^{a}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 1 
figure(5)
plot(tViL,Vpcca,'-','LineWidth', 1); hold all;
plot(tViL,iLa,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Load}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{Load}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 2 
figure(6)
plot(tViL,Vpcca,'-','LineWidth', 1); hold all;
plot(tViL,iLa,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Load}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{Load}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% iLb
figure(7)
plot(tViL,Vpccb,'-','LineWidth', 1); hold all;
plot(tViL,iLb,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Load}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{b}$','$i_{Load}^{b}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(8)
plot(tViL,Vpccb,'-','LineWidth', 1); hold all;
plot(tViL,iLb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Load}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{b}$','$i_{Load}^{b}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(9)
plot(tViL,Vpccb,'-','LineWidth', 1); hold all;
plot(tViL,iLb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Load}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQb = legend('$V_{pcc}^{b}$','$i_{Load}^{b}$');
% set(legend_PQb,'Interpreter','latex','Location','best');
hold off

% iLc
figure(10)
plot(tViL,Vpccc,'-','LineWidth', 1); hold all;
plot(tViL,iLc,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Load}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQc = legend('$V_{pcc}^{c}$','$i_{Load}^{c}$');
set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(11)
plot(tViL,Vpccc,'-','LineWidth', 1); hold all;
plot(tViL,iLc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Load}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{Load}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(12)
plot(tViL,Vpccc,'-','LineWidth', 1); hold all;
plot(tViL,iLc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Load}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{Load}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

%% Voltage grid currents comparison
clear all

load('TLSCViRedPhasesRobustPI_20200117.mat')

tViRed = TLSCVIPhasesRobustPI20200117(:,1);
Vpcca = TLSCVIPhasesRobustPI20200117(:,2);
iReda = TLSCVIPhasesRobustPI20200117(:,3);
Vpccb = TLSCVIPhasesRobustPI20200117(:,4);
iRedb = TLSCVIPhasesRobustPI20200117(:,5);
Vpccc = TLSCVIPhasesRobustPI20200117(:,6);
iRedc = TLSCVIPhasesRobustPI20200117(:,7);

% iReda
figure(13)
plot(tViRed,Vpcca,'-','LineWidth', 1); hold all;
plot(tViRed,iReda,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Red}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{a}$','$i_{Red}^{a}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 1 
figure(14)
plot(tViRed,Vpcca,'-','LineWidth', 1); hold all;
plot(tViRed,iReda,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Red}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{Red}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 2 
figure(15)
plot(tViRed,Vpcca,'-','LineWidth', 1); hold all;
plot(tViRed,iReda,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{Red}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{Red}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% iLb
figure(16)
plot(tViRed,Vpccb,'-','LineWidth', 1); hold all;
plot(tViRed,iRedb,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Red}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{b}$','$i_{Red}^{b}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(17)
plot(tViRed,Vpccb,'-','LineWidth', 1); hold all;
plot(tViRed,iRedb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Red}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{b}$','$i_{Red}^{b}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(18)
plot(tViRed,Vpccb,'-','LineWidth', 1); hold all;
plot(tViRed,iRedb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{Red}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQb = legend('$V_{pcc}^{b}$','$i_{Red}^{b}$');
% set(legend_PQb,'Interpreter','latex','Location','best');
hold off

% iLc
figure(19)
plot(tViRed,Vpccc,'-','LineWidth', 1); hold all;
plot(tViRed,iRedc,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Red}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQc = legend('$V_{pcc}^{c}$','$i_{Red}^{c}$');
set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(20)
plot(tViRed,Vpccc,'-','LineWidth', 1); hold all;
plot(tViRed,iRedc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Red}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{Red}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(21)
plot(tViRed,Vpccc,'-','LineWidth', 1); hold all;
plot(tViRed,iRedc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{Red}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{Red}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

%% Voltage SAPF currents comparison
clear all

load('TLSCVipPhasesRobustPI_20200117.mat')

tViS = TLSCVipPhasesRobustPI20200117(:,1);
Vpcca = TLSCVipPhasesRobustPI20200117(:,2);
iSa = TLSCVipPhasesRobustPI20200117(:,3);
Vpccb = TLSCVipPhasesRobustPI20200117(:,4);
iSb = TLSCVipPhasesRobustPI20200117(:,5);
Vpccc = TLSCVipPhasesRobustPI20200117(:,6);
iSc = TLSCVipPhasesRobustPI20200117(:,7);

% iReda
figure(22)
plot(tViS,Vpcca,'-','LineWidth', 1); hold all;
plot(tViS,iSa,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{S}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{a}$','$i_{S}^{a}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 1 
figure(23)
plot(tViS,Vpcca,'-','LineWidth', 1); hold all;
plot(tViS,iSa,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{S}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{S}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

%Zoom 2 
figure(24)
plot(tViS,Vpcca,'-','LineWidth', 1); hold all;
plot(tViS,iSa,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{a}$ v.s. $i_{S}^{a}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{a}$','$i_{S}^{a}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% iLb
figure(25)
plot(tViS,Vpccb,'-','LineWidth', 1); hold all;
plot(tViS,iSb,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{S}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQa = legend('$V_{pcc}^{b}$','$i_{S}^{b}$');
set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(26)
plot(tViS,Vpccb,'-','LineWidth', 1); hold all;
plot(tViS,iSb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{S}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.6 0.7])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQa = legend('$V_{pcc}^{b}$','$i_{S}^{b}$');
% set(legend_PQa,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(27)
plot(tViS,Vpccb,'-','LineWidth', 1); hold all;
plot(tViS,iSb,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{b}$ v.s. $i_{S}^{b}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQb = legend('$V_{pcc}^{b}$','$i_{S}^{b}$');
% set(legend_PQb,'Interpreter','latex','Location','best');
hold off

% iLc
figure(28)
plot(tViS,Vpccc,'-','LineWidth', 1); hold all;
plot(tViS,iSc,':','LineWidth', 2);
xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{S}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.5 2.5])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
legend_PQc = legend('$V_{pcc}^{c}$','$i_{S}^{c}$');
set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 1
figure(29)
plot(tViS,Vpccc,'-','LineWidth', 1); hold all;
plot(tViS,iSc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{S}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([0.8 0.9])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{S}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

% Zoom 2
figure(30)
plot(tViS,Vpccc,'-','LineWidth', 1); hold all;
plot(tViS,iSc,':','LineWidth', 2);
% xlabel('Time (seg)', 'Interpreter', 'Latex','FontSize',32);
% ylabel('$0.1*V_{pcc}^{c}$ v.s. $i_{S}^{c}$', 'Interpreter', 'Latex','FontSize',32);
xlim([2.0 2.1])
ylim([-20 20])
set(gca,'FontSize',32,'FontName','Times')
% legend_PQc = legend('$V_{pcc}^{c}$','$i_{S}^{c}$');
% set(legend_PQc,'Interpreter','latex','Location','best');
hold off

