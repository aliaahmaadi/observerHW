clc
clear
close all
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
 C = [1 1 1 1 1];
 D = 0;

%% reduced order observer 
A22 = [0 1 0 0;-16.2664 0 0 0.2531; 5.5413 0 0 -3.9805;0 0 -0.2255 -11.0345];
A12 = [0 0 1 0];
A21 = [0; 0; 0; 0];
A11 = 0;
B1 = 0;
B2 = [0 0 0 6.8966]';
C1 = 1;
C2 = [0 0 0 0];

dpoles = 9*[-2 -3 -4 -5];
L = place(A22',(C1*A12)',dpoles)';
F_o = A22 - L*C1*A12;
H_o = B2 - L*C1*B1;
G_o = inv(C1)*(A21 - L*C1*A11 + F_o*L*C1);



%% full state feedback controller
t = [-8 -1.23 -7 -3 -0.75];
 K = place(A,B,t);

k1 = K(1)
k2 = K(2)
k3 = K(3)
k4 = K(4)
k5 = K(5)
%% simulation
load_system('reduced_obs2')
       sim('reduced_obs2')
       
figure(1)
plot(ans.tout,ans.th','LineWidth',2)
hold on
plot(ans.tout,ans.th_o,'LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('\theta(radian)','LineWidth',2,'fontsize',14)
title('regulating \theta in presence of reduced order observer(\theta_0 = 20degree) ')
legend('plant','observer')
grid on
figure(2)
plot(ans.tout,ans.thd,'LineWidth',2)
hold on
plot(ans.tout,ans.thd_o,'LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('\rate of theta','LineWidth',2,'fontsize',14)
title('reduced order observer (rate of \theta)')
legend('plant','observer')
grid on
 figure(3)
plot(ans.tout,ans.phid,'LineWidth',2)
hold on
plot(ans.tout,ans.phid_o,'LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('rate of \phi','LineWidth',2,'fontsize',14)
title('reduced order observer')
legend('plant','observer')
grid on
figure(4)
plot(ans.tout,ans.i,'LineWidth',2)
hold on
plot(ans.tout,ans.i_o,'LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('current i','LineWidth',2,'fontsize',14)
title('reduced order observer')
legend('plant','observer') 
grid on
figure(5)
plot(ans.tout,ans.effort,'LineWidth',2)
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('amplitude (voltage)','LineWidth',2,'fontsize',14)
title('actuator effort')
grid on
figure(6)
plot(ans.tout,ans.phi,'LineWidth',2)


xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('\phi','LineWidth',2,'fontsize',14)
grid on