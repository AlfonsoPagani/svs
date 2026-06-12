close all
clear all
clc

%% INPUT
%
%                m4
%       ^        |
%     /PAY\      >
%    |_____|     m3
%    | ADP |     |
%    |_____|     >  k3
%    |     |     |
%    | FUE |     m2 
%    |_____|     |
%    |     |     <  k1
%    | O_X |     |
%    |_____|     m1
%      \ /
%      /_\
%
%       ^        F(t) = f, t>=0
%       |

m = [ 445e3
      116e3
      2e3
      15e3 ];                                                              % Discretized masses [kg]
  
k = [ 5e8
      5e8
      1e8 ];                                                               % Stiffness [N/m]
  
f = 7600e3 ;                                                               % Amplitude of step load [N] 775 ton

%% ANALYSIS

M_MAT = [ m(1) 0    0    0
          0    m(2) 0    0
          0    0    m(3) 0 
          0    0    0    m(4) ];                                           % Mass matrix
     
K_MAT = [ k(1) -k(1)       0           0
         -k(1)  k(1)+k(2) -k(2)        0
          0    -k(2)       k(2)+k(3)  -k(3) 
          0     0         -k(3)        k(3)];                              % Stiffness matrix
      
F_VEC = [ f
          0
          0
          0 ];
      
[PHI,lambda] = eig(K_MAT,M_MAT);                                           % Solve | K_MAT - omega^2 M_MAT | = 0 ... Normalize wrt mass


omega = real(sqrt(lambda));                                                % Natural periods [rad/s]
freq = omega ./ 2 ./ pi;                                                   % Natural frequencies [Hz]

M_GEN = PHI' * M_MAT * PHI;                                                % Generalized mass matrix

K_GEN = PHI' * K_MAT * PHI;                                                % Generalized stiffness matrix

F_GEN = PHI' * F_VEC;                                                      % Generalized force vector



t = linspace(0,160,1e5);                                                   % Discretized time vector

q1 = F_GEN(1) / 2 / M_GEN(1,1) * t.^2;                                     % Solution in terms of modal coordinates
q2 = F_GEN(2) / K_GEN(2,2) * ( 1 - cos(omega(2,2)*t) );
q3 = F_GEN(3) / K_GEN(3,3) * ( 1 - cos(omega(3,3)*t) );
q4 = F_GEN(4) / K_GEN(4,4) * ( 1 - cos(omega(4,4)*t) );

q = [ q1
      q2
      q3
      q4];
  
x = PHI * q;                                                               % Solution in terms of physical coordinates

q1dd = F_GEN(1) / M_GEN(1,1) .* ones(1,length(t));                         % Modal accelerations
q2dd = F_GEN(2) / K_GEN(2,2) * omega(2,2)^2 * cos(omega(2,2)*t);
q3dd = F_GEN(3) / K_GEN(3,3) * omega(3,3)^2 * cos(omega(3,3)*t);
q4dd = F_GEN(4) / K_GEN(4,4) * omega(4,4)^2 * cos(omega(4,4)*t);

x1dd = PHI(1,1) .* q1dd + PHI(1,2) .* q2dd + PHI(1,3) .* q3dd...
    + PHI(1,4) .* q4dd;    

x4dd = PHI(4,1) .* q1dd + PHI(4,2) .* q2dd + PHI(4,3) .* q3dd...
    + PHI(4,4) .* q4dd;                                                    % Payload acceleration

F4 = m(4) * x4dd;

%% POST-PROCESSING

for i=1:4
    PHIv(:,i)=PHI(:,i) ./ max(abs(PHI(:,i)));                              % normalize wrt max
end

pos = [ 0; 2; 4; 6];                                                       % Mode shapes
figure(1)
plot(zeros(4,1),pos,'k-o','LineWidth',3);
hold on
plot(zeros(4,1)+0.5,pos+PHIv(:,1),'r-o','LineWidth',3);
hold on
plot(zeros(4,1)+1,pos+PHIv(:,2),'r-o','LineWidth',3);
hold on
plot(zeros(4,1)+1.5,pos+PHIv(:,3),'r-o','LineWidth',3);
hold on
plot(zeros(4,1)+2.,pos+PHIv(:,4),'r-o','LineWidth',3);
xlim([-0.5 2]);
ylim('auto');
grid ON

figure(2)                                                                  % Lift
plot(t,x(1,:)./1e3,'-b','Linewidth',1.5)
xlabel('t [s]')
ylabel('Sollevamento [km]')
% saveas(gcf,'falconlift','epsc')

figure(3)                                                                  % Zoom lift
plot(t,x(1,:),'-b','Linewidth',1.5)
hold on
plot(t,x(2,:),'ob','Linewidth',1.,'MarkerIndices', 1:10:length(t))
hold on
plot(t,x(3,:),'--k','Linewidth',1.)
hold on 
plot(t,x(4,:),'sk','Linewidth',1.,'MarkerIndices', 1:50:length(t))
legend('Stg.1', 'Stg.2', 'Fair.+Adpt.', 'Payload','Location','northwest')
legend('boxoff')
xlabel('t [s]')
ylabel('Sollevamento [m]')
xlim([0 0.2])
% saveas(gcf,'falconliftzoom','epsc')

figure(4)                                                                  % Zoom payload acceleration
plot(t,x4dd./9.81,'-b','Linewidth',1)
xlabel('t [s]')
ylabel('Accelerazione [g]')
xlim([0 2])
yyaxis right
ylabel('Forza trasmessa [kN]')
xlim([0 2])
ylim([-4*9.81*m(4)/1e3 6*9.81*m(4)/1e3])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% saveas(gcf,'falconpayacc','epsc')

%% EFFECTIVE MODAL MASS

PHIR=[1;1;1;1];

LK = PHI' * M_MAT * PHIR;

for i=1:4 % Modal participation factor
    PARFACT(i) = LK(i) / M_GEN(i,i);
end

for i=1:4
   M_EFF(i) = LK(i)^2./M_GEN(i,i); % Modal effective mass
end

%% Free-Free Participation Factor

for i=2:4
    Dll(i-1)=trace((PHI(:,i)./lambda(i,i))*PHI(:,i)');
end

Dlltot=zeros(4,4);
for i=2:4
    Dlltot=Dlltot + (PHI(:,i)./lambda(i,i))*PHI(:,i)';
end
Dlltot=trace(Dlltot);

LFREE = Dll./Dlltot;

%%

x1m  = PHI(:,1) .* q1;
x12m = PHI(:,1) .* q1 + PHI(:,2) .* q2;
x123m = PHI(:,1) .* q1 + PHI(:,2) .* q2 + PHI(:,3) .* q3;

x4ddl = PHI(4,1) .* q1dd + PHI(4,2) .* q2dd;

figure(6)                                                                  % Zoom payload acceleration
plot(t,x4dd./9.81,'-b','Linewidth',1)
hold on
plot(t,x4ddl./9.81,'--r','Linewidth',1)
xlabel('t [s]')
ylabel('Accelerazione [g]')
% legend('Payload acc.')
xlim([0 2])
% ylim([-600 900])
ax = gca;
ax.YAxis(1).Color = 'k';
% saveas(gcf,'falconpayaccmodes12','epsc')

figure(7)                                                                  % Zoom lift
plot(t,x(4,:),'-b','Linewidth',1.5)
hold on
plot(t,x1m(4,:),'ob','Linewidth',1.,'MarkerIndices', 1:20:length(t))
hold on
plot(t,x12m(4,:),'--r','Linewidth',1.)
hold on
plot(t,x123m(4,:),'sk','Linewidth',1.,'MarkerIndices', 1:13:length(t))
legend('Exact', 'Mode 1', 'Modes 1 & 2', 'Modes 1, 2 & 3', 'Location','northwest')
legend('boxoff')
xlabel('t [s]')
ylabel('Sollevamento [m]')
xlim([0 0.2])
% saveas(gcf,'falconliftzoommodes','epsc')

