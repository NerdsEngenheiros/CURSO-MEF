clc
clear

#dados de entrada
L=3*10^3; #[mm]

E=200*10^3; %[MPa]

#A=(200*360)-(120*300);
A=(200*360)

#I=((200*360^3)/12)-((120*300^3)/12);
I=((200*360^3)/12)

 
#cargas externas
F4=(8*10^3)*cosd(45); #[N]
F5=-(8*10^3)*sind(45); #[N]
M6=0; #[Nmm]
Fr=[F4; F5; M6];

#matriz de rigidez
a=A*E/L;
b=E*I/(L^3);

k=[a 0 0 -a 0 0; 0 12*b 6*b*L 0 -12*b 6*b*L; 0 6*b*L 4*b*L^2 0 -6*b*L 2*b*L^2; -a 0 0 a 0 0; 0 -12*b -6*b*L 0 12*b -6*b*L; 0 6*b*L 2*b*L^2 0 -6*b*L 4*b*L^2];

#matriz de rigidez reduzida
kr=k(4:6,4:6);

#deslocamentos nodais
Dr=inv(kr)*Fr;
D=[0;0;0;Dr]

#forças nodais
F1=k(1,:)*D
F2=k(2,:)*D
M3=k(3,:)*D