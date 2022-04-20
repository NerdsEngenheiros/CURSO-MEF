clc
clear

E=200*10^3; %MPa = N/mm²

base=200; 
altura=200;

A=base*altura;
I=(base*altura^3)/12;

La=3*10^3;
Lb=1*10^3;
L=[La; Lb];


F4=20*10^3;
F5=0;
M6=0;
F7=0;
F8=-20*10^3;
M9=0;

Fr=[F4; F5; M6; F7; F8; M9];

a=(A*E)./L;
b=(E*I)./(L.^3);

alpha=[90; 0];
lamb=cosd(alpha);
mu=sind(alpha);
lamb2=lamb.^2;
mu2=mu.^2;
lambmu=lamb.*mu;

for n=1:length(L)
    K(:,:,n)=[a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n), -a(n)*lamb2(n)-12*b(n)*mu2(n),-(a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n); (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n), -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2; -a(n)*lamb2(n)-12*b(n)*mu2(n), -(a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n), a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n); -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n), (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2];
end

KG=zeros(9,9);

m=1;
for n=1:length(L)
KG(m:m+5,m:m+5)=KG(m:m+5,m:m+5)+K(:,:,n)
m=m+3    
end


%n=1
%Ka=[a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n), -a(n)*lamb2(n)-12*b(n)*mu2(n),-(a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n); (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n), -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2; -a(n)*lamb2(n)-12*b(n)*mu2(n), -(a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n), a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n); -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n), (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2];
%n=2
%Kb=[a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n), -a(n)*lamb2(n)-12*b(n)*mu2(n),-(a(n)-12*b(n))*lambmu(n), -6*b(n)*L(n)*mu(n); (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n), -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), 6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2; -a(n)*lamb2(n)-12*b(n)*mu2(n), -(a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n), a(n)*lamb2(n)+12*b(n)*mu2(n), (a(n)-12*b(n))*lambmu(n), 6*b(n)*L(n)*mu(n); -(a(n)-12*b(n))*lambmu(n), -a(n)*mu2(n)-12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n), (a(n)-12*b(n))*lambmu(n), a(n)*mu2(n)+12*b(n)*lamb2(n), -6*b(n)*L(n)*lamb(n); -6*b(n)*L(n)*mu(n), 6*b(n)*L(n)*lamb(n), 2*b(n)*L(n)^2, 6*b(n)*L(n)*mu(n), -6*b(n)*L(n)*lamb(n), 4*b(n)*L(n)^2];

%KG=zeros(9,9)
%KG(1:6,1:6)=Ka;
%KG(4:9,4:9)=KG(4:9,4:9)+Kb

Kr=KG(4:9,4:9)

Dr=inv(Kr)*Fr
D=[0;0;0;Dr]

F1=KG(1,:)*D
F2=KG(2,:)*D
M3=KG(3,:)*D
