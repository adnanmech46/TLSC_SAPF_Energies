% The following code computes the PI controllers based on the pole-placement method 
% presented in "Robust Control for Shunt Active Power Filters: 
% A Dynamical Model-based Approach with VerifiedControllability", 
% by Jorge-Humberto Urrea-Quintero, Jesús M. López-Lezama
% Nicolás Muñoz-Galeano, Energies, 2020.

% It is a .m file, and therefore it is meant to be used on Matlab.

% Code written by: Jorge-Humberto Urrea-Quintero.
% 2020

%%
clear all;
close all; clc;
%%
s=tf([1 0], [0 1]);
%System parameters

Vpcc = 120*sqrt(3);
f = 50;
w = 2*pi*f;

VL = 120*sqrt(3);
Vdc = 600;

fswp = 20e3;
Lp = 30.2e-3;
Rp = 1;

Cdc = 2200e-6;
Ro = 1000;

%% Non-linear model dq0

syms ipd ipq ip0 vdc Ev vpccd vpccq vpcc0 upd upq up0

dipd = (1/Lp)*(Lp*w*ipq-Rp*ipd-upd*(vdc/2)+vpccd);

dipq = (1/Lp)*(-Lp*w*ipd-Rp*ipq-upq*(vdc/2)+vpccq);

dip0 = (1/Lp)*(-Rp*ip0-up0*(vdc/2)-sqrt(3)*(Ev)+vpcc0);


dvdc = (1/Cdc)*(upd*ipd+upq*ipq+up0*ip0-(vdc/Ro));

dEv = (1/(2*Cdc))*(sqrt(3)*ip0-((2*Ev)/Ro));

%% Space-state representation

sistema = [dipd dipq dip0 dvdc dEv];
x = [ipd ipq ip0 vdc Ev];
u = [upd upq up0 vpccd vpccq vpcc0];
y = [ipd ipq ip0 vdc Ev];

Asyms = jacobian(sistema,x);       %matrix A simbolica
Bsyms = jacobian(sistema,u);       %matrix B simbolica
Csyms = jacobian(y,x);
Dsyms = jacobian(y,u);

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
%% Linearization

A = eval(Asyms);      %valor de la matriz A(sistema)
B = eval(Bsyms);      %valor de la matriz B
C = eval(Csyms);
D = eval(Dsyms);

sys = ss(A,B,C,D);

G = tf(sys);
Gidud = G(1,1);
Giduq = G(1,2);
Gidvpccd = G(1,4);
Gidvpccq = G(1,5);
Giquq = G(2,2);
Giqud = G(2,1);
Giqvpccd = G(2,4);
Giqvpccq = G(2,5);
Gi0u0 = G(3,3);
Gi0vpcc0 = G(3,6);


% sisotool(-Gidud)
% break

% % PI CONTROLLER
% Damping 0.707, Settling time = 0.3 ms
aid = 15452;
bid = 0.00011;
 
Tid = bid
Kid = aid*Tid
 
Gcid= -aid*(1+bid*s)/s; 

aiq = 14911;
biq = 0.00012;
 
Tiq = biq
Kiq = aiq*Tiq
 
Gciq= -aiq*(1+biq*s)/s; 

ai0 = 15122;
bi0 = 0.00012;

Ti0 = bi0
Ki0 = ai0*Ti0
 
Gci0= -ai0*(1+bi0*s)/s; 
%break

Sid = 1/(1+Gidud*Gcid);
Tidref = (Gidud*Gcid)/(1+Gidud*Gcid);
Tiduq = (Giduq)/(1+Gidud*Gcid);
Tidvpccd = (Gidvpccd)/(1+Gidud*Gcid);
Tidvpccq = (Gidvpccq)/(1+Gidud*Gcid);

Siq = 1/(1+Giquq*Gciq);
Tiqref = (Giquq*Gciq)/(1+Giquq*Gciq);
Tiqud = (Giqud)/(1+Giquq*Gciq);
Tiqvpccd = (Giqvpccd)/(1+Giquq*Gciq);
Tiqvpccq = (Giqvpccq)/(1+Giquq*Gciq);

Si0 = 1/(1+Gi0u0*Gci0);
Ti0ref = (Gi0u0*Gci0)/(1+Gi0u0*Gci0);
Ti0vpcc0 = (Gi0vpcc0)/(1+Gi0u0*Gci0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mag,phase,woutTidref] = bode(Tidref);
for i = 1:length(mag(1,1,:))
    magTidref(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTidref(i,1) = phase(1,1,i);
end

[mag,phase,woutTiduq] = bode(Tiduq);
for i = 1:length(mag(1,1,:))
    magTiduq(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTiduq(i,1) = phase(1,1,i);
end

[mag,phase,woutTidvpccd] = bode(Tidvpccd);
for i = 1:length(mag(1,1,:))
    magTidvpccd(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTidvpccd(i,1) = phase(1,1,i);
end

[mag,phase,woutTidvpccq] = bode(Tidvpccq);
for i = 1:length(mag(1,1,:))
    magTidvpccq(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTidvpccq(i,1) = phase(1,1,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mag,phase,woutTiqref] = bode(Tiqref);
for i = 1:length(mag(1,1,:))
    magTiqref(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTiqref(i,1) = phase(1,1,i);
end

[mag,phase,woutTiqud] = bode(Tiqud);
for i = 1:length(mag(1,1,:))
    magTiqud(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTiqud(i,1) = phase(1,1,i);
end

[mag,phase,woutTiqvpccd] = bode(Tiqvpccd);
for i = 1:length(mag(1,1,:))
    magTiqvpccd(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTiqvpccd(i,1) = phase(1,1,i);
end

[mag,phase,woutTiqvpccq] = bode(Tiqvpccq);
for i = 1:length(mag(1,1,:))
    magTiqvpccq(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTiqvpccq(i,1) = phase(1,1,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mag,phase,woutTi0ref] = bode(Ti0ref);
for i = 1:length(mag(1,1,:))
    magTi0ref(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTi0ref(i,1) = phase(1,1,i);
end

[mag,phase,woutTi0vpcc0] = bode(Ti0vpcc0);
for i = 1:length(mag(1,1,:))
    magTi0vpcc0(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTi0vpcc0(i,1) = phase(1,1,i);
end

% TiIo = (GiLIo)/(1+GiLd*Gci);

%% DC-link controller

Am = -1/(Ro*Cdc);
Bm = [upd/Cdc upq/Cdc];
Cm = 1;
Dm = [0 0];

SSvdc = ss(Am, Bm, Cm, Dm);
Gvdc = tf(SSvdc);

% sisotool(Gvdc)
% setting time = 0.02s
% damping factor = 0.704
% cut frequency = 30Hz

avdc = 29.441;
bvdc = 0.014;
 
Tvdc = bvdc
Kvdc = avdc*Tvdc
 
Gcvdc= avdc*(1+bvdc*s)/s;

Svdc = 1/(1+Gvdc(1,1)*Gcvdc);
Tvdcref = (Gvdc(1,1)*Gcvdc)/(1+Gvdc(1,1)*Gcvdc);
Tvdcuq = (Gvdc(1,2)*Gcvdc)/(1+Gvdc(1,2)*Gcvdc);

[mag,phase,woutTvdcref] = bode(Tvdcref);
for i = 1:length(mag(1,1,:))
    magTvdcref(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTvdcref(i,1) = phase(1,1,i);
end

[mag,phase,woutTvdcuq] = bode(Tvdcuq);
for i = 1:length(mag(1,1,:))
    magTvdcuq(i,1) = mag(1,1,i);
end
for i = 1:length(phase(1,1,:))
    phaseTvdcuq(i,1) = phase(1,1,i);
end

%% Sensitive functions
[mag,phase,woutSid] = bode(Sid);
for i = 1:length(mag(1,1,:))
    magSid(i,1) = mag(1,1,i);
end
Msid = max(magSid)

[mag,phase,woutSiq] = bode(Siq);
for i = 1:length(mag(1,1,:))
    magSiq(i,1) = mag(1,1,i);
end
Msiq = max(magSiq)

[mag,phase,woutSi0] = bode(Si0);
for i = 1:length(mag(1,1,:))
    magSi0(i,1) = mag(1,1,i);
end
Msi0 = max(magSi0)

[mag,phase,woutSvdc] = bode(Svdc);
for i = 1:length(mag(1,1,:))
    magSvdc(i,1) = mag(1,1,i);
end
Msid = max(magSvdc)
