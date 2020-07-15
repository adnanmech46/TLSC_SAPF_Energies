%% Calculo de las referencias para el SAPC en estado estacionario

clc; clear all; close all;
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