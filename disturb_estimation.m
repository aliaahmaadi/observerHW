clc
clear
%% System parameters
L=0.1810; % Artist's heigth
l=0.0476; % distance between artist's hands and wire
 % rod's lenght
M_M=0.51; % Artist's mass
M_R=0.39; % rod's mass
M_m=0.1;  %DC motor's %mass
M_e=0.075;  %encoder's mass
M_H=0;
J_M=0.0054; % artist's inertia
J_R=0.0488; % rod's inertia
J_H=0.1; % Housing's inertia
J_m=0.032; % motor's shaft inertia
j=9.7*10^-7; % gearbox's output shaft inertia
R_m=1.6; % DC motor's Electric resistance
L_m=0.145; % DC motor's inductance
K_m=0.0109; % motor's constant
N=3; % Gearbox transfer ratio
g=9.8; % gravity acceleration
%%


J=J_R+J_H+J_M;
J_RH=J_R+J_H;
M=M_m+M_R+M_H+M_e;

W=(J_M+M*l^2)*(J_RH+N^2*J_m+j)-(J_RH)^2;
Z=(J_M+M*l^2)*(J_RH+N^2*J_m+j)^2-J_RH^2*(J_RH+N^2*J_m+j);
T=(M_M*L/2+M*l)*(J_RH+N^2*J_m+j)*g;
H=-N*K_m*J_RH;
G=-J_RH*(M_M*L/2+M*l)*(J_RH+N^2*J_m+j)*g;
E=(J_RH^2+1)*(N*K_m);

A = [0 1 0 0 0
     T/W 0 0 0 H/W
     0 0 0 1 0
     G/Z 0 0 0 E/Z
     0 0 0 -(N*K_m)/L_m -R_m/L_m];
 B = [0 0 0 0 1/L_m]';
 C = [1 1 1 1 1]';
 D = [0 0 0 0 0]';
 
 
 %% pole placement
 t = [-8 -1.23 -7 -3 -0.75];
 K = place(A,B,t);

k1 = K(1)
k2 = K(2)
k3 = K(3)
k4 = K(4)
k5 = K(5)
kI = -150;
tt = 1.1*[-8 -1.23 -7 -3 -0.75];
 L = place(A,C,tt);

l1 = L(1)
l2 = L(2)
l3 = L(3)
l4 = L(4)
l5 = L(5)
l6 = 160;


%% simulation
load_system('disturbance_est')
       sim('disturbance_est')
       
figure(1)
plot(ans.tout,ans.th,'b','LineWidth',2)
%hold on
%plot(ans.tout,ans.th_h,'r','LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('\theta(radian)','LineWidth',2,'fontsize',14)
title('(Regulation \theta_0 = 20 degree) in presence of estimated disturbance')
%legend('plant','observer')
 
%figure(2)
%plot(ans.tout,ans.ef,'r','LineWidth',2)
%xlabel('Time','LineWidth',2,'fontsize',14)
%ylabel('amplitude (voltag)','LineWidth',2,'fontsize',14)
%title('actuator effort')
 
 
figure(2)
plot(ans.tout,ans.di,'b','LineWidth',2)
hold on
plot(ans.tout,ans.di_o,'r','LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('amplitude','LineWidth',2,'fontsize',14)
title('disturbance')
legend('plant','observer')
