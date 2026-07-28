clear all
close all
clc

%int diam 4, ext diam 6
l1=0.78035; 
l2=1.38325;
l3=1.14825;

alpha=deg2rad(22.5);

% step=[0.1 0.25 0.5 0.75 1 1.25 1.5 1.75 2];
step=[0.1 0.4 0.8 1.2 1.6 2];

for k=1:length(step)
    
    rA=step(k);    
    rB=rA*(cos(alpha)+sqrt((l1/rA)^2-(sin(alpha)^2)));
    rD=rA*(cos(alpha)+sqrt((l2/rA)^2-(sin(alpha)^2)));
    rC=rD*(cos(alpha)-sqrt((l3/rD)^2-(sin(alpha)^2)));
    
    figure(k)
    for i=1:1:8
        h=polarplot([2*alpha*(i-1) 2*alpha*(i-1)],[rA rC],'-b','LineWidth',3);
        hold on
        h=polarplot([alpha*(2*i-1) alpha*(2*i-1)],[rB rD],'-b','LineWidth',3);
        hold on
        h=polarplot([2*alpha*i 2*alpha*i],[rA rC],'-b','LineWidth',3);
        hold on
        h=polarplot([2*alpha*(i-1) alpha*(2*i-1)],[rA rB],'-ok','LineWidth',1);
        hold on
        h=polarplot([2*alpha*(i-1) alpha*(2*i-1)],[rA rD],'-ok','LineWidth',1);
        hold on
        h=polarplot([2*alpha*(i-1) alpha*(2*i-1)],[rC rD],'-ok','LineWidth',1);
        hold on
        h=polarplot([2*alpha*i alpha*(2*i-1)],[rA rB],'-ok','LineWidth',1);
        hold on
        h=polarplot([2*alpha*i alpha*(2*i-1)],[rA rD],'-ok','LineWidth',1);
        hold on
        h=polarplot([2*alpha*i alpha*(2*i-1)],[rC rD],'-ok','LineWidth',1);
        hold on
    end
    Ax = gca; % current axes
    Ax.ThetaAxisUnits = 'radians';
    thetaticks(0:pi/4:2*pi-pi/4)
    rlim([0 3])
    rticklabels({'r = 0','r = 2','r = 4','r = 6'})
    Ax.ThetaColor = 'blue';
    Ax.RColor = [0 .5 0];
    saveas(gcf,sprintf('DepSolArray_%i', k),'epsc')
    
    a1=rC-rA;
    a2=rD-rB;
    data(k,:)=[rA rB rC rD a1 a2];
    
end

figure(k+1)
plot(data(:,5),data(:,1))
hold on
plot(data(:,5),data(:,2))
hold on
plot(data(:,5),data(:,3))
hold on
plot(data(:,5),data(:,4))
legend('A','B','C','D')
xlabel('a1 [m]')
ylabel('Radial coordinate [m]')
grid on

figure(k+2)
plot(data(:,6),data(:,1))
hold on
plot(data(:,6),data(:,2))
hold on
plot(data(:,6),data(:,3))
hold on
plot(data(:,6),data(:,4))
legend('A','B','C','D')
xlabel('a2 [m]')
ylabel('Radial coordinate [m]')

