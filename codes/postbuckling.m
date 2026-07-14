clear all
close all
clc

%% Elastica
% Data una trave di alluminio, con modulo elastico E = 75 GPa, di
% lunghezza L = 0.5 m e sezione quadrata con lato b = 5 mm, calcolare: 
% * lo stato di equilibrio in post-buckling per 0<alpha<180deg
% * confrontare la soluzione nonlineare con la soluzione linearizzata

E = 75e9;
b = 0.005; J = b^4/12;
L = 0.5;

%% Equilibrio

alpha = linspace(0,pi,100);
rhosq = (sin(alpha/2)).^2;
K = ellipke(rhosq);
k = K./L;
P = k.^2*E*J;
y_L = 2*sqrt(rhosq)./k;

figure
plot(y_L,P)
title('Equilibrio')
xlabel('y(L), m')
ylabel('P, N')
hold on
plot([0 0.45], [38.55 38.55], 'r')
legend('post-buckling nl.', 'sol. linearizzata')

