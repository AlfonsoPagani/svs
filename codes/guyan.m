clear all
close all
clc

m1=445e3;
m2=116e3;
m3=3e3;
m4=4.5e3;
m5=3e3;
m6=1.5e3;
m7=1.5e3;
m8=1.5e3;
m9=1.5e3;
m10=1.5e3;

k1=5e8;
k2=5e8;
k3=1e8;
k4=1e8;
k5=.5e8;
k6=.5e8;
k7=1e8;
k8=1e8;
k9=.5e8;

M=diag([m1,m2,m3,m4,m5,m6,m7,m8,m9,m10]);

K=[k1 -k1 0 0 0 0 0 0 0 0;
    -k1 k1+k2 -k2 0 0 0 0 0 0 0;
    0 -k2 k2+k3 -k3 0 0 0 0 0 0;
    0 0 -k3 k3+k4+k7+k8 -k4 0 0 -k7 -k8 0;
    0 0 0 -k4 k4+k5 -k5 0 0 0 0
    0 0 0 0 -k5 k5 0 0 0 0;
    0 0 0 0 0 0 k6 -k6 0 0;
    0 0 0 -k7 0 0 -k6 k6+k7 0 0;
    0 0 0 -k8 0 0 0 0 k8+k9 -k9;
    0 0 0 0 0 0 0 0 -k9 k9];

% Mvinc=M; Kvinc=K;
M(:,1)=[]; M(1,:)=[];
K(:,1)=[]; K(1,:)=[];

% mm: 2,6,7,10
mm=[1,5,6,9]; % nodo 1 vincolato
ss=setdiff(1:length(K(:,1)),mm);

[PHI,lambda] = eig(K,M,'vector');
omega = real(sqrt(lambda));
freq = omega ./ 2 ./ pi;

M_GEN = PHI' * M * PHI;
K_GEN = PHI' * K * PHI;

G=-inv(K(ss,ss))*K(ss,mm);
T=[eye(length(mm)); G];

MGR=T'*[M(mm,mm),M(mm,ss);M(ss,mm),M(ss,ss)]*T;
KGR=T'*[K(mm,mm),K(mm,ss);K(ss,mm),K(ss,ss)]*T;

[PHI_m,lambda_GR] = eig(KGR,MGR,'vector');
[lambda_GR, ind] = sort(lambda_GR);
PHI_m = PHI_m(:, ind);
omega_GR = real(sqrt(lambda_GR));
freq_GR = omega_GR ./ 2 ./ pi;

PHI_GR=T*PHI_m;
PHI_GR_ok=zeros(length(K(:,1)),length(mm));
PHI_GR_ok(mm,:)=PHI_GR(1:length(mm),:);
PHI_GR_ok(ss,:)=PHI_GR(length(mm)+1:end,:);

M_GEN_GR = PHI_GR_ok' * M * PHI_GR_ok;
K_GEN_GR = PHI_GR_ok' * K * PHI_GR_ok;
for i=1:length(mm)
    PHI_GR_ok(:,i)=PHI_GR_ok(:,i)./M_GEN_GR(i,i);
end
M_GEN_GR = PHI_GR_ok' * M * PHI_GR_ok;
K_GEN_GR = PHI_GR_ok' * K * PHI_GR_ok;

N=length(mm);
MAC=zeros(N+1,N+1);
PHI_O=PHI(:,1:N);

for i=1:N
    for j=1:N
        NUM=PHI_O(:,i)'*PHI_GR_ok(:,j)*PHI_GR_ok(:,j)'*PHI_O(:,i);
        DEN=PHI_O(:,i)'*PHI_O(:,i)*PHI_GR_ok(:,j)'*PHI_GR_ok(:,j);
        MAC(i,j)=NUM/DEN;
    end
end

meanMAC=0;
for i=1:N
    meanMAC=meanMAC+MAC(i,i);
end
meanMAC=meanMAC/N
X=1:1:N;

figure(1);
pcolor(MAC);
Nc = 256; % number of colors
cmap = [linspace(1,0,Nc).' linspace(1,0,Nc).' ones(Nc,1)]; % decreasing R and G; B = 1
colormap(cmap)
t=colorbar;
axis square;
set(gca,'XTick',[1.5:1:N+0.5]);
set(gca,'XTickLabels',X);
xlabel('Modi di vibrare modello ORIGINALE'); %CUF
set(gca,'YTick',[1.5:1:N+0.5]);
set(gca,'YTickLabels',X);
ylabel('Modi di vibrare modello RIDOTTO'); %NAS
t.Limits = [0 1.0];
set(get(t,'ylabel'),'String', 'MAC');
print -depsc -r300 MAC2D.eps

CO=zeros(N+1,N+1);
for i=1:N
    for j=1:N
        NUM=PHI_O(:,i)'*M*PHI_GR_ok(:,j);
        DEN=PHI_O(:,i)'*M*PHI_O(:,i)*PHI_GR_ok(:,j)'*M*PHI_GR_ok(:,j);
        CO(i,j)=NUM/sqrt(DEN);
    end
end
CO=abs(CO);

figure(2);
pcolor(CO);
Nc = 256; % number of colors
cmap = [linspace(1,0,Nc).' linspace(1,0,Nc).' ones(Nc,1)]; % decreasing R and G; B = 1
colormap(cmap)
t=colorbar;
axis square;
set(gca,'XTick',[1.5:1:N+0.5]);
set(gca,'XTickLabels',X);
xlabel('Modi di vibrare modello ORIGINALE'); %CUF
set(gca,'YTick',[1.5:1:N+0.5]);
set(gca,'YTickLabels',X);
ylabel('Modi di vibrare modello RIDOTTO'); %NAS
t.Limits = [0 1.0];
set(get(t,'ylabel'),'String', 'COC');
print -depsc -r300 COC2D.eps

for i=1:length(PHI(1,:))
    PHImax(:,i)=PHI(:,i)./max(abs(PHI(:,i)));
end

for i=1:length(PHI_GR_ok(1,:))
    PHI_GR_okmax(:,i)=PHI_GR_ok(:,i)./max(abs(PHI_GR_ok(:,i)));
end
