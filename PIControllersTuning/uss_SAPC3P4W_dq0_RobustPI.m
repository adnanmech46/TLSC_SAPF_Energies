% The following code computes all the robust PI controllers 
% presented in "Robust Control for Shunt Active Power Filters: 
% A Dynamical Model-based Approach with VerifiedControllability", 
% by Jorge-Humberto Urrea-Quintero, Jesús M. López-Lezama
% Nicolás Muñoz-Galeano, Energies, 2020.

% It is a .m file, and therefore it is meant to be used on Matlab. The code
% employes the system control toolbox. Particularly, it relies on the systune
% function available in the lastest version of the software: v. R2020a.

% Code written by: Jorge-Humberto Urrea-Quintero.
% 2020

%% Uncertainty Model for a TLSC SAPF
clear all;
close all;
clc;
%% Uncertainty parameters

s=tf([1 0], [0 1]);
%Parametros del sistema

Vpcc = 120*sqrt(3);
f = 50;
w = 2*pi*f;

VL = 120*sqrt(3);
Vdc = 600;

fswp = 20e3;

% Lp = 30.2; % ---> to ms
% Rp = 1;
% 
% Cdc = 2200e-3; % ---> to ms
% Ro = 1000;

Lp = ureal('Lp',30.2e-3,'Percentage',[-30, 30]);
Rp = ureal('Rp',1,'Percentage',[-50, 50]);

Cdc = ureal('Cdc',2200e-6,'Percentage',[-30, 30]);
Ro = ureal('Ro',1000,'Percentage',[-50, 50]);

% Lp = ureal('Lp',30.2e-3,'Percentage',[-90, 90]);
% Rp = ureal('Rp',1,'Percentage',[-90, 90]);
% 
% Cdc = ureal('Cdc',2200e-6,'Percentage',[-90, 90]);
% Ro = ureal('Ro',1000,'Percentage',[-90, 90]);

%% Operating point

vpccd = sqrt(2)*120;
vpccq = 0;
vpcc0 = 0;

ipd = 1.18163308543578;
ipq = 4.37418877873929;
ip0 = 0;

vdc = 600;
Ev = 0;

upd = 0.700081969280322;
upq = -0.0519502081215973;
up0 = 0;

%% State-space model

A = ...
    [  -Rp/Lp,  100*pi,               0, -upd/(2*Lp),           0;...
    -100*pi,  -Rp/Lp,               0, -upq/(2*Lp),           0;...
    0,       0,          -Rp/Lp, -up0/(2*Lp), -3^(1/2)/Lp;...
    upd/Cdc, upq/Cdc,         up0/Cdc, -1/(Cdc*Ro),           0;...
    0,       0, 3^(1/2)/(2*Cdc),           0, -1/(Cdc*Ro)];

B = ...
    [ -(2500*vdc)/151,               0,               0, 5000/151,        0,        0;...
    0, -(2500*vdc)/151,               0,        0, 5000/151,        0;...
    0,               0, -(2500*vdc)/151,        0,        0, 5000/151;...
    (5000*ipd)/11,   (5000*ipq)/11,   (5000*ip0)/11,        0,        0,        0;...
    0,               0,               0,        0,        0,        0];

C = ...
    [ 1, 0, 0, 0, 0;...
    0, 1, 0, 0, 0;...
    0, 0, 1, 0, 0;...
    0, 0, 0, 1, 0;...
    0, 0, 0, 0, 1];

D = ...
    [ 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0];

uss1 = ss(A,B,C,D);
% G = tf(uss1);
% [Gnom1,DeltaNom1,Blkstruct,Normunc] = lftdata(uss1);
SSnom1 = getNominal(uss1);

%% PI Tuning
% R1 = TuningGoal.Sensitivity('dLoad',tf([1.20 0],[1 4e3]));
% % viewGoal(R1)
% R2 = TuningGoal.MaxLoopGain('dLoad',10,4e3);
% % viewGoal(R2)
% R3 = TuningGoal.Poles('dLoad',0,0.707,25e3);
% % viewGoal(R3)
% R4 = TuningGoal.Margins('dLoad',20,60);
% % viewGoal(R4)

R1 = TuningGoal.Sensitivity('dLoad',tf([1.20 0],[1 4e3]));
% viewGoal(R1)
R2 = TuningGoal.MaxLoopGain('dLoad',10,4e3);
% viewGoal(R2)
R3 = TuningGoal.Poles('dLoad',0,0.707,3*4e3);
% viewGoal(R3)
R4 = TuningGoal.Margins('dLoad',15,60);
% viewGoal(R4)

SoftReqs = [R1 R2 R3];
HardReqs = [R4];

% SoftReqs = [R1 R2];
% HardReqs = [R3 R4];

% id

CPIid = tunablePID('C','pi');
APPIid = AnalysisPoint('dLoad');
CL0PIid = feedback(APPIid*SSnom1(1,1)*CPIid,1);
CL0PIid.InputName = 'idRef';
CL0PIid.OutputName = 'ud';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIid,fSoftPIid] = systune(CL0PIid,SoftReqs,HardReqs,opt);
showTunable(CLPIid)

SPIid = getSensitivity(CLPIid,'dLoad');
figure
step(getNominal(SPIid),3e-3)
title('Load disturbance rejection')
legend('Nominal')

figure
viewGoal(R1,CLPIid)
figure
viewGoal(R2,CLPIid)
figure
viewGoal(R3,CLPIid)
figure
viewGoal(R4,CLPIid)

% iq

CPIiq = tunablePID('C','pi');
APPIiq = AnalysisPoint('dLoad');
CL0PIiq = feedback(APPIiq*SSnom1(2,2)*CPIiq,1);
CL0PIiq.InputName = 'iqRef';
CL0PIiq.OutputName = 'uq';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIiq,fSoftPIiq] = systune(CL0PIiq,SoftReqs,HardReqs,opt);
showTunable(CLPIiq)

SPIiq = getSensitivity(CLPIiq,'dLoad');
figure
step(getNominal(SPIiq),3e-3)
title('Load disturbance rejection')
legend('Nominal')

figure
viewGoal(R1,CLPIid)
figure
viewGoal(R2,CLPIid)
figure
viewGoal(R3,CLPIiq)
figure
viewGoal(R4,CLPIiq)

% i0

CPIi0 = tunablePID('C','pi');
APPIi0 = AnalysisPoint('dLoad');
CL0PIi0 = feedback(APPIi0*SSnom1(3,3)*CPIi0,1);
CL0PIi0.InputName = 'i0Ref';
CL0PIi0.OutputName = 'u0';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIi0,fSoftPIi0] = systune(CL0PIi0,SoftReqs,HardReqs,opt);
showTunable(CLPIi0)

SPIi0 = getSensitivity(CLPIi0,'dLoad');
figure
step(getNominal(SPIi0),3e-3)
title('Load disturbance rejection')
legend('Nominal')

figure
viewGoal(R1,CLPIi0)
figure
viewGoal(R2,CLPIi0)
figure
viewGoal(R3,CLPIi0)
figure
viewGoal(R4,CLPIi0)

%% Robust PI Tuning
% id

CPIidRobust = tunablePID('C','pi');
APPIidRobust = AnalysisPoint('dLoad');
CL0PIidRobust = feedback(APPIidRobust*uss1(1,1)*CPIidRobust,1);
CL0PIidRobust.InputName = 'idRef';
CL0PIidRobust.OutputName = 'ud';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIidRobust,fSoftPIidRobust] = systune(CL0PIidRobust,SoftReqs,HardReqs,opt);
showTunable(CLPIidRobust)

SPIidRobust = getSensitivity(CLPIidRobust,'dLoad');
figure
step(usample(SPIidRobust,30),getNominal(SPIidRobust),3e-3)
title('Load disturbance rejection')
legend('Sampled uncertainty','Nominal')

figure
viewGoal(R1,CLPIidRobust)
figure
viewGoal(R2,CLPIidRobust)
figure
viewGoal(R3,CLPIidRobust)
figure
viewGoal(R4,CLPIidRobust)

% iq

CPIiqRobust = tunablePID('C','pi');
APPIiqRobust = AnalysisPoint('dLoad');
CL0PIiqRobust = feedback(APPIiqRobust*uss1(2,2)*CPIiqRobust,1);
CL0PIiqRobust.InputName = 'iqRef';
CL0PIiqRobust.OutputName = 'uq';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIiqRobust,fSoftPIiqRobust] = systune(CL0PIiqRobust,SoftReqs,HardReqs,opt);
showTunable(CLPIiqRobust)

SPIiqRobust = getSensitivity(CLPIiqRobust,'dLoad');
figure
step(usample(SPIiqRobust,30),getNominal(SPIiqRobust),3e-3)
title('Load disturbance rejection')
legend('Sampled uncertainty','Nominal')

figure
viewGoal(R1,CLPIiqRobust)
figure
viewGoal(R2,CLPIiqRobust)
figure
viewGoal(R3,CLPIiqRobust)
figure
viewGoal(R4,CLPIiqRobust)

% i0

CPIi0Robust = tunablePID('C','pi');
APPIi0Robust = AnalysisPoint('dLoad');
CL0PIi0Robust = feedback(APPIi0Robust*uss1(3,3)*CPIi0Robust,1);
CL0PIi0Robust.InputName = 'i0Ref';
CL0PIi0Robust.OutputName = 'u0';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIi0Robust,fSoftPIi0Robust] = systune(CL0PIi0Robust,SoftReqs,HardReqs,opt);
showTunable(CLPIi0Robust)

SPIi0Robust = getSensitivity(CLPIi0Robust,'dLoad');
figure
step(usample(SPIi0Robust,30),getNominal(SPIi0Robust),3e-3)
title('Load disturbance rejection')
legend('Sampled uncertainty','Nominal')

figure
viewGoal(R1,CLPIi0Robust)
figure
viewGoal(R2,CLPIi0Robust)
figure
viewGoal(R3,CLPIi0Robust)
figure
viewGoal(R4,CLPIi0Robust)

%% Transfer functions

% L = Lp.^(1/2);
% Ipq = ipq;
% RL = Rp;
% 
% Ipd = (1/(2*RL))*(vpccd-sqrt((((vpccd^2)-4*RL*((RL*(Ipq^2))+((vdc^2)/(2*Ro)))))));
% % Ipd = (1/(2*RL))*(vpccd-(((vpccd^2)-4*RL*((RL*(Ipq^2))+((vdc^2)/(2*Ro)))))^(0.5));
% Upd = (1/Ipd)*((vdc/Ro)+(2*Ipq/vdc)*(p.L*w*Ipd+RL*Ipq));
% Upq = -(2/vdc)*(p.L*w*Ipd+RL*Ipq);
% 
% GidudMaple = (-2*Cdc*L*Ro*Vdc*s^2+(-2*Cdc*Ro*Rp*Vdc-2*Ipd*L*Ro*Upd-2*L*Vdc)*s-2*Ipd*L*Ro*Upq*w-2*Ipd*Ro*Rp*Upd-Ro*Upq^2*Vdc-2*Rp*Vdc)...
%     /(4*Cdc*L^2*Ro*s^3+(8*Cdc*L*Ro*Rp+4*L^2)*s^2+(4*Cdc*L^2*Ro*w^2+4*Cdc*Ro*Rp^2+2*L*Ro*Upd^2+2*L*Ro*Upq^2+8*L*Rp)*s+4*w^2*L^2+2*Ro*Rp*Upd^2+2*Ro*Rp*Upq^2+4*Rp^2);
% 
% GiduqMaple = ((-2*Cdc*L*Ro*Vdc*w-2*Ipq*L*Ro*Upd)*s-2*Ipq*L*Ro*Upq*w-2*Ipq*Ro*Rp*Upd+Ro*Upd*Upq*Vdc-2*L*w*Vdc)...
%     /(4*Cdc*L^2*Ro*s^3+(8*Cdc*L*Ro*Rp+4*L^2)*s^2+(4*Cdc*L^2*Ro*w^2+4*Cdc*Ro*Rp^2+2*L*Ro*Upd^2+2*L*Ro*Upq^2+8*L*Rp)*s+4*w^2*L^2+2*Ro*Rp*Upd^2+2*Ro*Rp*Upq^2+4*Rp^2);
% 
% GiqudMaple = ((2*Cdc*L*Ro*Vdc*w-2*Ipd*L*Ro*Upq)*s+2*Ipd*L*Ro*Upd*w-2*Ipd*Ro*Rp*Upq+Ro*Upd*Upq*Vdc+2*L*w*Vdc)...
%     /(4*Cdc*L^2*Ro*s^3+(8*Cdc*L*Ro*Rp+4*L^2)*s^2+(4*Cdc*L^2*Ro*w^2+4*Cdc*Ro*Rp^2+2*L*Ro*Upd^2+2*L*Ro*Upq^2+8*L*Rp)*s+4*w^2*L^2+2*Ro*Rp*Upd^2+2*Ro*Rp*Upq^2+4*Rp^2);
% 
% GiquqMaple = (-2*Cdc*L*Ro*Vdc*s^2+(-2*Cdc*Ro*Rp*Vdc-2*Ipq*L*Ro*Upq-2*L*Vdc)*s+2*Ipq*L*Ro*Upd*w-2*Ipq*Ro*Rp*Upq-Ro*Upd^2*Vdc-2*Rp*Vdc)...
%     /(4*Cdc*L^2*Ro*s^3+(8*Cdc*L*Ro*Rp+4*L^2)*s^2+(4*Cdc*L^2*Ro*w^2+4*Cdc*Ro*Rp^2+2*L*Ro*Upd^2+2*L*Ro*Upq^2+8*L*Rp)*s+4*w^2*L^2+2*Ro*Rp*Upd^2+2*Ro*Rp*Upq^2+4*Rp^2);
% 
% GvdcudMaple = (2*Ipd*L^2*Ro*s^2+(4*Ipd*L*Rp-L*Upd*Vdc)*Ro*s+(2*Ipd*L^2*w^2+L*Upq*Vdc*w+2*Ipd*Rp^2-Rp*Upd*Vdc)*Ro)...
%     /(2*Cdc*L^2*Ro*s^3+(4*Cdc*L*Ro*Rp+2*L^2)*s^2+(2*Cdc*L^2*Ro*w^2+2*Cdc*Ro*Rp^2+L*Ro*Upd^2+L*Ro*Upq^2+4*L*Rp)*s+2*w^2*L^2+Ro*Rp*Upd^2+Ro*Rp*Upq^2+2*Rp^2);
% 
% GvdcuqMaple = (2*Ipq*L^2*Ro*s^2+(4*Ipq*L*Rp-L*Upq*Vdc)*Ro*s+(2*Ipq*L^2*w^2-L*Upd*Vdc*w+2*Ipq*Rp^2-Rp*Upq*Vdc)*Ro)...
%     /(2*Cdc*L^2*Ro*s^3+(4*Cdc*L*Ro*Rp+2*L^2)*s^2+(2*Cdc*L^2*Ro*w^2+2*Cdc*Ro*Rp^2+L*Ro*Upd^2+L*Ro*Upq^2+4*L*Rp)*s+2*w^2*L^2+Ro*Rp*Upd^2+Ro*Rp*Upq^2+2*Rp^2);
% 
% % Zero secuence
% 
% Gi0u0Maple = (-Cdc*Ro*Vdc*s-Vdc)/(2*Cdc*L*Ro*s^2+(2*Cdc*Ro*Rp+2*L)*s+3*Ro+2*Rp);

%% Decomposing uncertain object

% [Gnom1,DeltaNom1,Blkstruct,Normunc] = lftdata(uss1);

%% Input uncertainty
% W1 = makeweight(0.1,6300,10);
% W2 = makeweight(0.1,6300,10);
% W3 = makeweight(0.1,6300,10);
% W4 = makeweight(0.25,630,10);
% W5 = makeweight(0.25,630,10);
% W6 = makeweight(0.25,630,10);
% Delta1 = ultidyn('Delta1',[1 1]);
% Delta2 = ultidyn('Delta2',[1 1]);
% Delta3 = ultidyn('Delta2',[1 1]);
% Delta4 = ultidyn('Delta2',[1 1]);
% Delta5 = ultidyn('Delta2',[1 1]);
% Delta6 = ultidyn('Delta2',[1 1]);
%
% %
% % Uncertainty model
% W = blkdiag(W1,W2,W3,W4,W5,W6);
% figure(1)
% omega = logspace(0,6,2000);
% bodemag(W1,'r-',W2,'b--',W3,'g',omega), grid
% title('Magnitude responses of weighting filters')
% legend('W_1','W_2','W_3')
% Delta = blkdiag(Delta1,Delta2,Delta3,Delta4,Delta5,Delta6);
% G = uss1*(eye(6) + Delta*W);

%% Controller transfer matrix
% id

% PI CONTROLLER
% Damping 0.707, Settling time = 0.3 ms
aid = 15452;
bid = 0.00011;

Tid = bid
Kid = aid*Tid

Gcid= -aid*(1+bid*s)/s;

CPIidPP = tunablePID('C','pi');
CPIidPP.Kp.Value = Kid;        
CPIidPP.Ki.Value = Kid/Tid; 
CPIidPP.Kp.Free = false;
CPIidPP.Ki.Free = false;
APPIidPP = AnalysisPoint('dLoad');
CL0PIidPP = feedback(APPIidPP*uss1(1,1)*CPIidPP,1,+1);
CL0PIidPP.InputName = 'idRef';
CL0PIidPP.OutputName = 'ud';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIidPP,fSoftPIidPP] = systune(CL0PIidPP,SoftReqs,HardReqs,opt);
showTunable(CLPIidPP)

SPIidPP = getSensitivity(CLPIidPP,'dLoad');

APPIidUncertain = AnalysisPoint('dLoad');
% CLPIidNominal = -0.631-(1.19e+03/s);
CLPIidNominal = CLPIid.Blocks.C.Kp.Value + CLPIid.Blocks.C.Ki.Value/s;
CLPIidUncertain = feedback(APPIidUncertain*uss1(1,1)*CLPIidNominal,1);
SPIidUncertain = getSensitivity(CLPIidUncertain,'dLoad');

figure
step(getNominal(SPIid),usample(SPIidUncertain,30),getNominal(SPIidRobust),usample(SPIidRobust,30),getNominal(SPIidPP),usample(SPIidPP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
step(getNominal(CLPIid),usample(CLPIidUncertain,30),getNominal(CLPIidRobust),usample(CLPIidRobust,30),getNominal(CLPIidPP),usample(CLPIidPP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
viewGoal(R1,CLPIidPP)
figure
viewGoal(R2,CLPIidPP)
figure
viewGoal(R3,CLPIidPP)
figure
viewGoal(R4,CLPIidPP)

% iq

aiq = 14911;
biq = 0.00012;

Tiq = biq
Kiq = aiq*Tiq

Gciq= -aiq*(1+biq*s)/s;

CPIiqPP = tunablePID('C','pi');
CPIiqPP.Kp.Value = Kiq;        
CPIiqPP.Ki.Value = Kiq/Tiq; 
CPIiqPP.Kp.Free = false;
CPIiqPP.Ki.Free = false;
APPIiqPP = AnalysisPoint('dLoad');
CL0PIiqPP = feedback(APPIiqPP*uss1(2,2)*CPIiqPP,1,+1);
CL0PIiqPP.InputName = 'iqRef';
CL0PIiqPP.OutputName = 'uq';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIiqPP,fSoftPIiqPP] = systune(CL0PIiqPP,SoftReqs,HardReqs,opt);
showTunable(CLPIiqPP)

SPIiqPP = getSensitivity(CLPIiqPP,'dLoad');

APPIiqUncertain = AnalysisPoint('dLoad');
% CLPIiqNominal = -6.74-(0.252/s);
CLPIiqNominal = CLPIiq.Blocks.C.Kp.Value + CLPIiq.Blocks.C.Ki.Value/s;
CLPIiqUncertain = feedback(APPIiqUncertain*uss1(2,2)*CLPIiqNominal,1);
SPIiqUncertain = getSensitivity(CLPIiqUncertain,'dLoad');

figure
step(getNominal(SPIiq),usample(SPIiqUncertain,30),getNominal(SPIiqRobust),usample(SPIiqRobust,30),getNominal(SPIiqPP),usample(SPIiqPP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
step(getNominal(CLPIiq),usample(CLPIiqUncertain,30),getNominal(CLPIiqRobust),usample(CLPIiqRobust,30),getNominal(CLPIiqPP),usample(CLPIiqPP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
viewGoal(R1,CLPIiqPP)
figure
viewGoal(R2,CLPIiqPP)
figure
viewGoal(R3,CLPIiqPP)
figure
viewGoal(R4,CLPIiqPP)


% i0

ai0 = 15122;
bi0 = 0.00012;

Ti0 = bi0
Ki0 = ai0*Ti0

Gci0= -ai0*(1+bi0*s)/s;

CPIi0PP = tunablePID('C','pi');
CPIi0PP.Kp.Value = Ki0;        
CPIi0PP.Ki.Value = Ki0/Ti0; 
CPIi0PP.Kp.Free = false;
CPIi0PP.Ki.Free = false;
APPIi0PP = AnalysisPoint('dLoad');
CL0PIi0PP = feedback(APPIi0PP*uss1(3,3)*CPIi0PP,1,+1);
CL0PIi0PP.InputName = 'i0Ref';
CL0PIi0PP.OutputName = 'u0';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIi0PP,fSoftPIi0PP] = systune(CL0PIi0PP,SoftReqs,HardReqs,opt);
showTunable(CLPIi0PP)

SPIi0PP = getSensitivity(CLPIi0PP,'dLoad');

APPIi0Uncertain = AnalysisPoint('dLoad');
CLPIi0Nominal = CLPIi0.Blocks.C.Kp.Value + CLPIi0.Blocks.C.Ki.Value/s;
CLPIi0Uncertain = feedback(APPIi0Uncertain*uss1(3,3)*CLPIi0Nominal,1);
SPIi0Uncertain = getSensitivity(CLPIi0Uncertain,'dLoad');

figure
step(getNominal(SPIi0),usample(SPIi0Uncertain,30),getNominal(SPIi0Robust),usample(SPIi0Robust,30),getNominal(SPIi0PP),usample(SPIi0PP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
step(getNominal(CLPIi0),usample(CLPIi0Uncertain,30),getNominal(CLPIi0Robust),usample(CLPIi0Robust,30),getNominal(CLPIi0PP),usample(CLPIi0PP,30),3e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
viewGoal(R1,CLPIi0PP)
figure
viewGoal(R2,CLPIi0PP)
figure
viewGoal(R3,CLPIi0PP)
figure
viewGoal(R4,CLPIi0PP)

%% Performance comparison

TCLPIidPP = getIOTransfer(CLPIidPP,'idRef','ud');
TCLPIid = getIOTransfer(CLPIid,'idRef','ud');
TCLPIidRobust = getIOTransfer(CLPIidRobust,'idRef','ud');

figure
clf
subplot(131), wcsigma(TCLPIidPP,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('PP tuning')
subplot(132), wcsigma(TCLPIid,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Nominal tuning'), legend('off')
subplot(133), wcsigma(TCLPIidRobust,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Robust tuning'), legend('off')

TCLPIiqPP = getIOTransfer(CLPIiqPP,'iqRef','uq');
TCLPIiq = getIOTransfer(CLPIiq,'iqRef','uq');
TCLPIiqRobust = getIOTransfer(CLPIiqRobust,'iqRef','uq');

figure
clf
subplot(131), wcsigma(TCLPIiqPP,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('PP tuning')
subplot(132), wcsigma(TCLPIiq,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Nominal tuning'), legend('off')
subplot(133), wcsigma(TCLPIiqRobust,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Robust tuning'), legend('off')

TCLPIi0PP = getIOTransfer(CLPIi0PP,'i0Ref','u0');
TCLPIi0 = getIOTransfer(CLPIi0,'i0Ref','u0');
TCLPIi0Robust = getIOTransfer(CLPIi0Robust,'i0Ref','u0');

figure
clf
subplot(131), wcsigma(TCLPIi0PP,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('PP tuning')
subplot(132), wcsigma(TCLPIi0,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Nominal tuning'), legend('off')
subplot(133), wcsigma(TCLPIi0Robust,{1e1,1e6}), grid
set(gca,'YLim',[-10 2.5]), title('Robust tuning'), legend('off')

%% MIMO controller

% k11 = Gcid;
% k22 = Gciq;
% k33 = Gci0;

% K  = blkdiag(k11,k22,k33);
% 
% feedin = [1 2 3];
% feedout = [1 2 3];
% CLPIidPolePlacement = feedback(SSnom1,K,feedin,feedout);
% CL0PIidPolePlacement = feedback(SSnom1*K,1);

% Determine sensitivity functions
% looptransfer = loopsens(uss1(1:3,1:3),K);
% Ti = looptransfer.Ti;
% To = looptransfer.To;

%% Robust stability analysis with robuststab
% omega = logspace(-1,6,200);
% To_g = ufrd(To,omega);
% opt = robopt('Display','on');
% [stabmarg,destabunc,report,info] = robuststab(To_g,opt);
% stabmarg
% 
% report
% 
% pole(usubs(To,destabunc))
% 
% figure(2)
% semilogx(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
% grid
% title('Robust stability')
% xlabel('Frequency (rad/s)')
% ylabel('mu')
% legend('\mu-upper bound','\mu-lower bound','NorthWest')

%% DC-link loop

Am = -1/(Ro*Cdc);
Bm = [upd/Cdc upq/Cdc];
Cm = 1;
Dm = [0 0];

SSvdc = ss(Am, Bm, Cm, Dm);
Gvdc = tf(SSvdc);

% Opt PI

R1vdc = TuningGoal.Sensitivity('dLoad',tf([1.20 0],[1 2e3/5]));
% viewGoal(R1vdc)
R2vdc = TuningGoal.MaxLoopGain('dLoad',10,8e3/5);
% viewGoal(R2vdc)
R3vdc = TuningGoal.Poles('dLoad',0,0.707,25e3/250);
% viewGoal(R3vdc)
R4vdc = TuningGoal.Margins('dLoad',20,60);
% viewGoal(R4vdc)

% SoftReqsvdc = [R1vdc R2vdc R3vdc];
% HardReqsvdc = [R4vdc];

SoftReqsvdc = [R1vdc R2vdc];
HardReqsvdc = [R3vdc R4vdc];

CPIvdc = tunablePID('C','pi');
APPIvdc = AnalysisPoint('dLoad');
CL0PIvdc = feedback(APPIvdc*getNominal(SSvdc(1,1))*CPIvdc,1);
CL0PIvdc.InputName = 'iLRef';
CL0PIvdc.OutputName = 'iL';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIvdc,fSoftPIvdc] = systune(CL0PIvdc,SoftReqsvdc,HardReqsvdc,opt);
showTunable(CLPIvdc)

SPIvdc = getSensitivity(CLPIvdc,'dLoad');
figure
step(getNominal(SPIvdc),150e-3)
title('Load disturbance rejection')
legend('Nominal')

figure
viewGoal(R1vdc,CLPIvdc)
figure
viewGoal(R2vdc,CLPIvdc)
figure
viewGoal(R3vdc,CLPIvdc)
figure
viewGoal(R4vdc,CLPIvdc)

% Robust PI

CPIvdcRobust = tunablePID('C','pi');
APPIvdcRobust = AnalysisPoint('dLoad');
CL0PIvdcRobust = feedback(APPIvdcRobust*SSvdc(1,1)*CPIvdcRobust,1);
CL0PIvdcRobust.InputName = 'iLRef';
CL0PIvdcRobust.OutputName = 'iL';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIvdcRobust,fSoftPIvdcRobust] = systune(CL0PIvdcRobust,SoftReqsvdc,HardReqsvdc,opt);
showTunable(CLPIvdcRobust)

SPIvdcRobust = getSensitivity(CLPIvdcRobust,'dLoad');
figure
step(usample(SPIvdcRobust,30),getNominal(SPIvdcRobust),150e-3)
title('Load disturbance rejection')
legend('Sampled uncertainty','Nominal')

figure
viewGoal(R1vdc,CLPIvdcRobust)
figure
viewGoal(R2vdc,CLPIvdcRobust)
figure
viewGoal(R3vdc,CLPIvdcRobust)
figure
viewGoal(R4vdc,CLPIvdcRobust)

% PP - PI

% sisotool(Gvdc)
% setting time = 0.02s
% damping factor = 0.704
% cut frequency = 30Hz

avdc = 29.441;
bvdc = 0.014;

Tvdc = bvdc
Kvdc = avdc*Tvdc

Gcvdc= avdc*(1+bvdc*s)/s;

CPIvdcPP = tunablePID('C','pi');
CPIvdcPP.Kp.Value = Kvdc;        
CPIvdcPP.Ki.Value = Kvdc/Tvdc; 
CPIvdcPP.Kp.Free = false;
CPIvdcPP.Ki.Free = false;
APPIvdcPP = AnalysisPoint('dLoad');
CL0PIvdcPP = feedback(APPIvdcPP*SSvdc(1,1)*CPIvdcPP,1);
CL0PIvdcPP.InputName = 'iLRef';
CL0PIvdcPP.OutputName = 'iL';

opt = systuneOptions('RandomStart',2);
rng(0);
[CLPIvdcPP,fSoftPIvdcPP] = systune(CL0PIvdcPP,SoftReqsvdc,HardReqsvdc,opt);
showTunable(CLPIvdcPP)

SPIvdcPP = getSensitivity(CLPIvdcPP,'dLoad');
figure
step(usample(SPIvdcPP,30),getNominal(SPIvdcPP),150e-3)
title('Load disturbance rejection')
legend('Sampled uncertainty','Nominal')

APPIvdcUncertain = AnalysisPoint('dLoad');
CLPIvdcNominal = CLPIvdc.Blocks.C.Kp.Value + CLPIvdc.Blocks.C.Ki.Value/s;
CLPIvdcUncertain = feedback(APPIvdcUncertain*SSvdc(1,1)*CLPIvdcNominal,1);
SPIvdcUncertain = getSensitivity(CLPIvdcUncertain,'dLoad');

figure
step(getNominal(SPIvdc),usample(SPIvdcUncertain,30),getNominal(SPIvdcRobust),usample(SPIvdcRobust,30),getNominal(SPIvdcPP),usample(SPIvdcPP,30),150e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
step(getNominal(CLPIvdc),usample(CLPIvdcUncertain,30),getNominal(CLPIvdcRobust),usample(CLPIvdcRobust,30),getNominal(CLPIvdcPP),usample(CLPIvdcPP,30),150e-3)
title('Load disturbance rejection')
legend('Nominal - OptPI', 'Uncertain - OptPI','Nominal - RobustPI', 'Uncertain - RobustPI','Nominal - PI-PP', 'Uncertain - PI-PP')

figure
viewGoal(R1vdc,CLPIvdcPP)
figure
viewGoal(R2vdc,CLPIvdcPP)
figure
viewGoal(R3vdc,CLPIvdcPP)
figure
viewGoal(R4vdc,CLPIvdcPP)

%% Performance comparison

TCLPIvdcPP = getIOTransfer(CLPIvdcPP,'iLRef','iL');
TCLPIvdc = getIOTransfer(CLPIvdc,'iLRef','iL');
TCLPIvdcRobust = getIOTransfer(CLPIvdcRobust,'iLRef','iL');

figure
clf
subplot(131), wcsigma(TCLPIvdcPP,{1e-2,1e3}), grid
set(gca,'YLim',[-10 3]), title('PP tuning')
subplot(132), wcsigma(TCLPIvdc,{1e-2,1e3}), grid
set(gca,'YLim',[-10 3]), title('Nominal tuning'), legend('off')
subplot(133), wcsigma(TCLPIvdcRobust,{1e-2,1e3}), grid
set(gca,'YLim',[-10 3]), title('Robust tuning'), legend('off')

% Determine sensitivity functions
% looptransferDC = loopsens(SSvdc(1,1),Gcvdc);
% TiDC = looptransferDC.Ti;
% ToDC = looptransferDC.To;

%% Robust stability analysis with robuststab
% omega = logspace(-1,6,200);
% To_gDC = ufrd(ToDC,omega);
% opt = robopt('Display','on');
% [stabmargDC,destabuncDC,reportDC,infoDC] = robuststab(To_gDC,opt);
% stabmargDC
% 
% reportDC
% 
% pole(usubs(ToDC,destabuncDC))
% 
% figure(3)
% semilogx(infoDC.MussvBnds(1,1),'r-',infoDC.MussvBnds(1,2),'b--')
% grid
% title('Robust stability')
% xlabel('Frequency (rad/s)')
% ylabel('mu')
% legend('\mu-upper bound','\mu-lower bound','NorthWest')
