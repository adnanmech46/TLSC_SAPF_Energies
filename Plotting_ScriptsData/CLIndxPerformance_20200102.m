% The following code computes all the performance indices 
% reported in "Robust Control for Shunt Active Power Filters: 
% A Dynamical Model-based Approach with VerifiedControllability", 
% by Jorge-Humberto Urrea-Quintero, Jesús M. López-Lezama
% Nicolás Muñoz-Galeano, Energies, 2020.

% It is a .m file, and therefore it is meant to be used on Matlab.
% To compute the indices, it is required to import the data on the .txt
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

idIAEPIPP = trapz(tPIPP,abs(idPIPPRef - idPIPP))
iqIAEPIPP = trapz(tPIPP,abs(iqPIPPRef - iqPIPP))
i0IAEPIPP = trapz(tPIPP,abs(i0PIPPRef - i0PIPP))
vdcIAEPIPP = trapz(tPIPP,abs(vdcPIPPRef - vdcPIPP))

idITAEPIPP = trapz(tPIPP,tPIPP.*abs(idPIPPRef - idPIPP))
iqITAEPIPP = trapz(tPIPP,tPIPP.*abs(iqPIPPRef - iqPIPP))
i0ITAEPIPP = trapz(tPIPP,tPIPP.*abs(i0PIPPRef - i0PIPP))
vdcITAEPIPP = trapz(tPIPP,tPIPP.*abs(vdcPIPPRef - vdcPIPP))

% PI - Opt

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

idIAEPIOpt = trapz(tPIOpt,abs(idPIOptRef - idPIOpt))
iqIAEPIOpt = trapz(tPIOpt,abs(iqPIOptRef - iqPIOpt))
i0IAEPIOpt = trapz(tPIOpt,abs(i0PIOptRef - i0PIOpt))
vdcIAEPIOpt = trapz(tPIOpt,abs(vdcPIOptRef - vdcPIOpt))

idITAEPIOpt = trapz(tPIOpt,tPIOpt.*abs(idPIOptRef - idPIOpt))
iqITAEPIOpt = trapz(tPIOpt,tPIOpt.*abs(iqPIOptRef - iqPIOpt))
i0ITAEPIOpt = trapz(tPIOpt,tPIOpt.*abs(i0PIOptRef - i0PIOpt))
vdcITAEPIOpt = trapz(tPIOpt,tPIOpt.*abs(vdcPIOptRef - vdcPIOpt))

% PI - Robust

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

idIAEPIRobust = trapz(tPIRobust,abs(idPIRobustRef - idPIRobust))
iqIAEPIRobust = trapz(tPIRobust,abs(iqPIRobustRef - iqPIRobust))
i0IAEPIRobust = trapz(tPIRobust,abs(i0PIRobustRef - i0PIRobust))
vdcIAEPIRobust = trapz(tPIRobust,abs(vdcPIRobustRef - vdcPIRobust))

idITAEPIRobust = trapz(tPIRobust,tPIRobust.*abs(idPIRobustRef - idPIRobust))
iqITAEPIRobust = trapz(tPIRobust,tPIRobust.*abs(iqPIRobustRef - iqPIRobust))
i0ITAEPIRobust = trapz(tPIRobust,tPIRobust.*abs(i0PIRobustRef - i0PIRobust))
vdcITAEPIRobust = trapz(tPIRobust,tPIRobust.*abs(vdcPIRobustRef - vdcPIRobust))

%% Control effort

upd = 0.700081969280322;
upq = -0.0519502081215973;
up0 = 0;
idref = 1.18163308543578;

% PI - Pole Placement

load('PIPP_ControlAction_20200102.mat')

tCAPP = PIPPControlAction20200102(:,1);

u0PP = PIPPControlAction20200102(:,2);
udPP = PIPPControlAction20200102(:,3);
uqPP = PIPPControlAction20200102(:,4);


IAUu0PP = trapz(tCAPP,abs(up0-u0PP))
IAUudPP = trapz(tCAPP,abs(upd-udPP))
IAUuqPP = trapz(tCAPP,abs(upq-uqPP))
IAUidrefPP = trapz(tCAPP,abs(idref-idPIPPRef))

% PI - Opt

load('PIOpt_ControlAction_20200102.mat')

tCAOpt = PIOptControlAction20200102(:,1);

u0Opt = PIOptControlAction20200102(:,2);
udOpt = PIOptControlAction20200102(:,3);
uqOpt = PIOptControlAction20200102(:,4);

IAUu0Opt = trapz(tCAOpt,abs(up0-u0Opt))
IAUudOpt = trapz(tCAOpt,abs(upd-udOpt))
IAUuqOpt = trapz(tCAOpt,abs(upq-uqOpt))
IAUidrefOpt = trapz(tCAOpt,abs(idref-idPIOptRef))

% PI - Robust

load('PIRobust_ControlAction_20200102.mat')

tCARobust = PIRobustControlAction20200102(:,1);

u0Robust = PIRobustControlAction20200102(:,2);
udRobust = PIRobustControlAction20200102(:,3);
uqRobust = PIRobustControlAction20200102(:,4);

IAUu0Robust = trapz(tCARobust,abs(up0-u0Robust))
IAUudRobust = trapz(tCARobust,abs(upd-udRobust))
IAUuqRobust = trapz(tCARobust,abs(upq-uqRobust))
IAUidrefRobust = trapz(tCARobust,abs(idref-idPIRobustRef))
