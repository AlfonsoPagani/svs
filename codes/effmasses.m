%********************************************************************
%******               MUL2 - STRUCTURAL DYNAMICS               ******
%******     Modal effective mass and partecipation factors     ******
%********************************************************************

clear all
close all
clc

%% DEFINIZIONE DEL PROBLEMA
%Definiamo il sistema in termini di matrici di rigidezza e massa;
m1 = 445e3; m2 = 116e3; m3 = 3e3;   m4 = 4.5e3; m5 = 3e3;
m6 = 1.5e3; m7 = 1.5e3; m8 = 1.5e3; m9=  1.5e3; m10= 1.5e3;

k1=5e8;  k2=5e8; k3=1e8; k4=1e8; k5=.5e8;
k6=.5e8; k7=1e8; k8=1e8; k9=.5e8;

M=diag([m1,m2,m3,m4,m5,m6,m7,m8,m9,m10]);
 
K=[k1   -k1     0       0           0      0     0    0     0     0;
  -k1  k1+k2   -k2      0           0      0     0    0     0     0;
    0   -k2   k2+k3    -k3          0      0     0    0     0     0;
    0    0     -k3  k3+k4+k7+k8   -k4      0     0   -k7   -k8    0;
    0    0      0      -k4       k4+k5   -k5     0    0     0     0
    0    0      0        0        -k5     k5     0    0     0     0;
    0    0      0        0          0      0    k6   -k6    0     0;
    0    0      0      -k7          0      0   -k6  k6+k7   0     0;
    0    0      0      -k8          0      0     0    0   k8+k9 -k9;
    0    0      0        0          0      0     0    0    -k9   k9];

DOF = length(K(1,:));

% Applico la condizione al contorno per i modi elastici. Imporre uno
% spostamento nullo significa rimuovere riga e colonna delle matrici di
% massa e rigidezza, per cui bloccando x1 considero la matrice [2:10,2:10]

K_BC = K(2:10,2:10);
M_BC = M(2:10,2:10);

ActiveDOF = length(K_BC(1,:));

%% ANALISI MODALE
% Problema agli autovalori: definisco e risolvo (K-w^2M)Phi = 0
[Phi, Lambda] = eig(K_BC,M_BC);

% Dagli autovalori lambda ricavo le frequenze naturali: fn = 1/2pi rad(lambda)
Frequencies = diag(sqrt(Lambda))/(2*pi);

fprintf('<strong>Mode #</strong>    <strong>Frequency [Hz]</strong>\n');
fprintf('-------------------------------------\n')
for i = 1:length(Frequencies)
    fprintf('% 5d      %8.6f\n', i,Frequencies(i));
end

% Osserviamo la soluzione in termini di autovettori ottenuti
fprintf('\n\n')
fprintf('<strong>Eigenvectors</strong>\n');
fprintf('----------------------------------------------------------\n')

% Autovettori riscalati rispetto al proprio massimo
for i = 1:size(Phi,1)
    Phi(:,i) = Phi(:,i)/max(Phi(:,i));
end

for i = 1:size(Phi,1)
    fprintf('%8.4f', Phi(i,:));
    fprintf('\n');
end

%% PLOT DEGLI AUTOVETTORI
% Plottiamo un generico autovettore per osservare la deformata
Psi = zeros(DOF,DOF);
Psi(2:DOF,2:DOF) = Phi;

ModeToPlot1 = 1;
ModeToPlot2 = 2;
ModeToPlot3 = 3;

Phi_Plot = [Psi(:,ModeToPlot1+1), Psi(:,ModeToPlot2+1), Psi(:,ModeToPlot3+1)];
Omega = [Frequencies(ModeToPlot1), Frequencies(ModeToPlot2),Frequencies(ModeToPlot3)];   % 1x3

animate_mode_shape(Phi_Plot, 0.7*Omega, 'Periods', 6, 'FPS', 120, ...
                   'FigurePosition', [100 200 1000 500]);



%% MASSE MODALI EFFETTIVE, FATTORI DI PARTECIPAZIONE

% I modi elastici sono i secondi 9, scriviamo quindi:
Phi_Rigido = diag(ones(ActiveDOF));

% Costruiamo le matrici generalizzate:
GeneralizedMass = (Phi') * M_BC * Phi;
GeneralizedStiffness = (Phi') * K_BC * Phi;

% Ripulisco valori spuri della matrice tendenti a zero macchina
for i=1:ActiveDOF
    for j = 1:ActiveDOF
        if(abs(GeneralizedMass(i,j))<1e-5)
            GeneralizedMass(i,j)=0;
        end
        if(abs(GeneralizedStiffness(i,j))<1e-5)
            GeneralizedStiffness(i,j)=0;
        end
    end
end

% Guardiamo le matrici
fprintf('\n\n')
fprintf('<strong>Generalized stiffness matrix</strong>\n');
fprintf('----------------------------------------------------------\n')

for i = 1:size(GeneralizedStiffness,1)
    fprintf('%12.3e', GeneralizedStiffness(i,:));
    fprintf('\n');
end

fprintf('\n\n')
fprintf('<strong>Generalized mass matrix</strong>\n');
fprintf('----------------------------------------------------------\n')

for i = 1:size(GeneralizedMass,1)
    fprintf('%12.3f', GeneralizedMass(i,:));
    fprintf('\n');
end


% Fattore di partecipazione modale:

L = Phi' * M_BC * Phi_Rigido;

ModalEffectiveMass = diag(zeros(ActiveDOF));
ModalParticipationFactor = diag(zeros(ActiveDOF));

for i = 1:ActiveDOF
    ModalParticipationFactor(i) = L(i) / GeneralizedMass(i,i);
    ModalEffectiveMass(i)       = L(i)^2 / GeneralizedMass(i,i);
end

fprintf('\n\n')
fprintf(['<strong>DOF</strong>     <strong>Partecipation factor  ' ...
    '</strong>     <strong>Modal Effective Mass</strong>\n']);
fprintf('----------------------------------------------------------\n')

for i = 1:ActiveDOF
    fprintf('%3d  %20.4f  %20.4f', i,ModalParticipationFactor(i),ModalEffectiveMass(i));
    fprintf('\n');
end

fprintf('----------------------------------------------------------\n')
fprintf('                        Total mass: %8.3f', sum(ModalEffectiveMass(1:DOF-1)));




