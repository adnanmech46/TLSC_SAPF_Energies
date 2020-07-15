% The following code computes the values of the states at the operating point 
% reported in "Robust Control for Shunt Active Power Filters: 
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

RL = 1;
f = 50;
w = 2*pi*f;
vdc = 600;
Ro = 1000;
L = 30.2e-3;
vpd = sqrt(2)*120;
ipq = -4.37418877873929

ipd = (1/(2*RL))*(vpd-sqrt((((vpd^2)-4*RL*((RL*(ipq^2))+((vdc^2)/(2*Ro)))))))
upd = (1/ipd)*((vdc/Ro)+(2*ipq/vdc)*(L*w*ipd+RL*ipq))
upq = -(2/vdc)*(L*w*ipd+RL*ipq)