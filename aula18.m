clc
clear

%DADOS DE ENTRADA
E=200*10^9; % [Pa]

L1=0.5; % [m]
L2=1.5; % [m]
L=[L1; L2];

b=0.2; % [m]
h=0.05; % [m]

I=(b*h^3)/12; 

%VETOR DE CARGAS NODAIS REDUZIDO
M2=0;
F3=-1*10^3;
M4=0;
M6=0;

FR=[M2; F3; M4; M6];

%MATRIZ DE RIGIDEZ
n=1;
k1=[12*E*I/L(n)^3 6*E*I/L(n)^2 -12*E*I/L(n)^3 6*E*I/L(n)^2; 6*E*I/L(n)^2 4*E*I/L(n) -6*E*I/L(n)^2 2*E*I/L(n); -12*E*I/L(n)^3 -6*E*I/L(n)^2 12*E*I/L(n)^3 -6*E*I/L(n)^2; 6*E*I/L(n)^2 2*E*I/L(n) -6*E*I/L(n)^2 4*E*I/L(n)];

n=2;
k2=[12*E*I/L(n)^3 6*E*I/L(n)^2 -12*E*I/L(n)^3 6*E*I/L(n)^2; 6*E*I/L(n)^2 4*E*I/L(n) -6*E*I/L(n)^2 2*E*I/L(n); -12*E*I/L(n)^3 -6*E*I/L(n)^2 12*E*I/L(n)^3 -6*E*I/L(n)^2; 6*E*I/L(n)^2 2*E*I/L(n) -6*E*I/L(n)^2 4*E*I/L(n)];

%MATRIZ DE RIGIDEZ GLOBAL
K=zeros(6,6);

K(1:4,1:4)=k1;

K(3:6,3:6)=K(3:6,3:6)+k2;

%MATRIZ DE RIGIDEZ REDUZIDO
KR=K([2:4,6],[2:4,6]);

%VETOR DE DESLOCAMENTO
DR=inv(KR)*FR;

D=[0; DR(1:3); 0 ; DR(4)]

%FORÇAS NODAIS
F1=K(1,:)*D
F5=K(5,:)*D


