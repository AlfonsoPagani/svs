clear all
close all
clc

%% Buckling Euleriano
% Data una trave di alluminio, con modulo elastico E = 75 GPa, di
% lunghezza L = 0.5 m e sezione quadrata con lato b = 5 mm, calcolare: 
% * i primi tre valori di carico critico
% * i primi tre modi di buckling

E = 75e9;
b = 0.005; J = b^4/12;
L = 0.5;

%% Carico critico

for i=1:3
    P(i) = (2*i-1)^2 * E * J * pi^2 / L^2 / 4;
end

figure
for i=1:3
    line([0 1],[P(i) P(i)])
    hold on
end
title('Equilibrio linearizzato')
xlabel('y, m')
ylabel('P, N)')

%% Modi di buckling
% y = B sin(kx)

x = 0:0.01:0.5;
y = zeros(length(x),3);
B = 1;
for i=1:3
    k = sqrt(P(i)/E/J);
    y(:,i) = B * sin(k*x);
end

figure
for i=1:3
    plot(x,y(:,i))
    hold on
end
legend('Modo 1','Modo 2','Modo 3')
title('Modi di buckling')
xlabel('x, m')
ylabel('y, m')