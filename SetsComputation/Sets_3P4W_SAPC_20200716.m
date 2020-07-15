% The following code computes the sets required for the controllability
% analysis reported in "Robust Control for Shunt Active Power Filters: 
% A Dynamical Model-based Approach with VerifiedControllability", 
% by Jorge-Humberto Urrea-Quintero, Jesús M. López-Lezama
% Nicolás Muñoz-Galeano, Energies, 2020.

% It is a .m file, and therefore it is meant to be used on Matlab.

% Code written by: Jorge-Humberto Urrea-Quintero.
% 2020

%%
clear all;
close all;
clc;
%% 
%Sytem parameters

Vpcc = 120*sqrt(3);
f = 50;
w = 2*pi*f;

VL = 120*sqrt(3);
Vdc = 600;
Ra = 18.58;
Rb = 18.58;
Rc = 18.58;
La = 44e-3;
Lb = 44e-3;
Lc = 44e-3;

fswp = 20e3;
Lp = 30.2e-3;
Rp = 1;

Cdc = 2200e-6;
Ro = 1000;

%Operating point

ipd = 1.18163308543578;
ipq = 4.37418877873929;
ip0 = 0;

vdc = 600;
Ev = 0;

upd = 0.700081969280322;
upq = -0.0519502081215973;
up0 = 0;

vpccd = sqrt(2)*120;
vpccq = 0;
vpcc0 = 0;

% Max and min states values
x_min = [-10;-10;550];
x_max = [10;10;650];

% States initial condition
x0=[ipd ipq vdc];

%Simulation time
t_ini = 0; 
t_fin = 50e-3;
Nm = 100000; % Sample size

%Control actions
U1min=0.1*upd;
U1max=2*upd;
U2min=-upq*0.1;
U2max=-upq*2;

%Non robust case
% U3min=sqrt(2)*120;
% U3max=sqrt(2)*120;

%Robust case
U3min=0.8*sqrt(2)*120;
U3max=1.2*sqrt(2)*120;

U1 = unifrnd(U1min,U1max,Nm,1);
U2 = -unifrnd(U2min,U2max,Nm,1);
U3 = unifrnd(U3min,U3max,Nm,1);

%Variables allocation
X_sol_p=x0;
X_sol_n=x0;

%% Positivo

% parpool(2)

tic
    parfor j=1:Nm
        U=[U1(j) U2(j) U3(j)];
        [t_traj x_traj]=ode45(@dynamicalSAPC3P4W,[t_ini:1/1000:t_fin],x0,[],U);
        for k=1:length(t_traj)
            if x_traj(k,1) > x_max(1) || x_traj(k,1) < x_min(1) || ...
                x_traj(k,2) > x_max(2) || x_traj(k,2) < x_min(2) ||...
                x_traj(k,3) > x_max(3) || x_traj(k,3) < x_min(3)
            x_traj(:,1)=x0(1)*ones(size(x_traj(:,1)));
            x_traj(:,2)=x0(2)*ones(size(x_traj(:,2)));
            x_traj(:,3)=x0(3)*ones(size(x_traj(:,3)));
            end
        end
        X_sol_p=[X_sol_p;x_traj(end,:)];
    end
toc
X_sol_p(isnan(X_sol_p(:,1)),1)=x0(1);
X_sol_p(isnan(X_sol_p(:,2)),2)=x0(2);
X_sol_p(isnan(X_sol_p(:,3)),2)=x0(3);

%% Negativo
tic
    parfor j=1:Nm
        U=[U1(j) U2(j) U3(j)];
        [t_traj x_traj]=ode45(@dynamicalSAPC3P4W,[t_fin:-1/1000:t_ini],x0,[],U);
        for k=1:length(t_traj)
            if x_traj(k,1) > x_max(1) || x_traj(k,1) < x_min(1) || ...
                x_traj(k,2) > x_max(2) || x_traj(k,2) < x_min(2) ||...
                x_traj(k,3) > x_max(3) || x_traj(k,3) < x_min(3)
%             x_traj(:,1)=x0(1)*ones(size(x_traj(:,1)));
%             x_traj(:,2)=x0(2)*ones(size(x_traj(:,2)));
%             x_traj(:,3)=x0(3)*ones(size(x_traj(:,3)));
            x_traj(:,1)=x0(1)
            x_traj(:,2)=x0(2)
            x_traj(:,3)=x0(3)
            end
        end
        X_sol_n=[X_sol_n;x_traj];
    end
toc
X_sol_n(isnan(X_sol_n(:,1)),1)=x0(1);
X_sol_n(isnan(X_sol_n(:,2)),2)=x0(2);
X_sol_n(isnan(X_sol_n(:,3)),2)=x0(3);

delete(gcp)

%% Calculo area
X_ro_p=[round(10.*X_sol_p(:,1)),round(10.*X_sol_p(:,2)),round(10.*X_sol_p(:,3))];     %Redondea Ca a 3 cifras, y T a 3 cifras
X_ro_n=[round(10.*X_sol_n(:,1)),round(10.*X_sol_n(:,2)),round(10.*X_sol_n(:,3))];     %Redondea Ca a 3 cifras, y T a 3 cifras

nrever=0;
nalcan=0;
ncontro=0;
tic
for i=0:10000                                                                 %Recorre todas las temperaturas
    pos_p=find(X_ro_p(:,3)==i);                                                 %Encuentra la posicion en las que está la temperatura actual
    pos_n=find(X_ro_n(:,3)==i);                                                 %Encuentra la posicion en las que está la temperatura actual
    X_rever=intersect(X_ro_p(pos_p,1:2),X_ro_n(pos_n,1:2));                         %Encuentra la intersección entre Ca y el vector de comparacion
    X_alcan=intersect(X_ro_p(pos_p,1:2),X_ro_p(pos_p,1:2));
    X_contro=intersect(X_ro_n(pos_n,1:2),X_ro_n(pos_n,1:2));
    nrever=nrever+numel(X_rever);                                               %Cuenta y almacena la cantidad de veces que coincidio
    nalcan=nalcan+numel(X_alcan);
    ncontro=ncontro+numel(X_contro);
end
toc

Acontro=ncontro
Aalcan=nalcan
Arever=nrever

% save resultnomd

% save resultoptd

% save resultpeord

%% 

figure(1)
plot(X_sol_n(:,2),X_sol_n(:,1),'.r')
hold on
plot(X_sol_p(:,2),X_sol_p(:,1),'.')
plot(x0(2),x0(1),'sk')
xlabel('i_{q}')
ylabel('i_{d}')
legend('Controllable set', 'Reachable set', 'Operating pint')

figure(2)
plot(X_sol_n(:,1),X_sol_n(:,3),'.r')
hold on
plot(X_sol_p(:,1),X_sol_p(:,3),'.')
plot(x0(1),x0(3),'sk')
xlabel('i_{d}')
ylabel('v_{DC}')
legend('Controllable set', 'Reachable set', 'Operating pint')

figure(3)
plot(X_sol_n(:,2),X_sol_n(:,3),'.r')
hold on
plot(X_sol_p(:,2),X_sol_p(:,3),'.')
plot(x0(2),x0(3),'sk')
xlabel('i_{q}')
ylabel('v_{DC}')
legend('Controllable set', 'Reachable set', 'Operating pint')

figure(4)
plot3(X_sol_n(:,1),X_sol_n(:,2),X_sol_n(:,3),'.r')
hold on
plot3(X_sol_p(:,1),X_sol_p(:,2),X_sol_p(:,3),'.')
plot3(x0(1),x0(2),x0(3),'sk')
grid on
xlabel('i_{d}')
ylabel('i_{q}')
zlabel('v_{DC}')
legend('Controllable set', 'Reachable set', 'Operating pint')