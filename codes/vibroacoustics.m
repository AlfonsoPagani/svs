clear all
close all
clc

M=0.5;
g=9.81;
z=0.05;
k=1.5*10^5;
fn=sqrt(k/M)/2/pi;
A=1;

pref=2*10^(-5);

SPL=[31.5   112
    63     123
    125    126
    250    135
    500    138
    1000   127
    2000   120];

df=(2^(1/2)-2^(-1/2)).*SPL(:,1);
f = SPL(:,1);

Wp=(pref*10.^(SPL(:,2)/20)).^2./(df(:));

H=abs(1./(1-2.*1i.*z.*(fn./f)-(fn./f).^2))*A./M;

Wa=H.^2.*Wp./(g.^2);

figure(1)
stairs(f-df./2,Wa,'LineWidth',2); grid on
set(gca,'XScale','log');
xlabel('f [Hz]'), ylabel('W_a [g^2/Hz]');

acc_rms = 0;
for i=1:length(Wa)
   acc_rms = acc_rms + Wa(i)*df(i); 
end
acc_rms = sqrt(acc_rms)

Wd = Wa.*(g.^2)./((2*pi*f).^4);

figure(2)
stairs(f-df./2,Wd,'LineWidth',2); grid on
set(gca,'XScale','log');
xlabel('f [Hz]'), ylabel('W_d [m^2/Hz]');

disp_rms = 0;
for i=1:length(Wd)
   disp_rms = disp_rms + Wd(i)*df(i); 
end
disp_rms = sqrt(disp_rms)