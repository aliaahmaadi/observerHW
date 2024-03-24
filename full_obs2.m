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
kI = 20;
tt = 1.1*[-8 -1.23 -7 -3 -0.75];
 L = place(A,C,tt);

l1 = L(1)
l2 = L(2)
l3 = L(3)
l4 = L(4)
l5 = L(5)
l6 = 100;
%% LQR method
%R = 100;
%Q = 0.01*eye(5,5);
%[ko,mo,e] = lqr(A,B,Q,R);

%k1 = ko(1)
%k2 = ko(2)
%k3 = ko(3)
%k4 = ko(4)
%k5 = ko(5)

%ei = eig(A-B*ko)
%tt = 3*[ei(1) ei(2) ei(3) ei(4) ei(5)];
% L = place(A,C,tt);

%l1 = L(1)
%l2 = L(2)
%l3 = L(3)
%l4 = L(4)
%l5 = L(5)

%% simulation

load_system('full_observer')
       sim('full_observer')
       
figure(1)
plot(ans.tout,ans.thh,'b')
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('\theta(radian)','LineWidth',2,'fontsize',14)
title('(Regulation \theta_0 = 20 degree) in presence of measurement noise')

 
figure(2)
plot(ans.tout,ans.eff,'k')
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('amplitude (voltag)','LineWidth',2,'fontsize',14)
title('actuator effort')
 
 
figure(3)
plot(ans.tout,ans.noi,'r')
xlabel('Time','LineWidth',2,'fontsize',14)
ylabel('amplitude','LineWidth',2,'fontsize',14)
title('measured noise')

% 
%figure(4)
%plot(ans.tout,ans.i_p,'b','LineWidth',2)
%hold on
%plot(ans.tout,ans.i_ob,'r','LineWidth',2)
%xlabel('Time','LineWidth',2,'fontsize',14)
%ylabel('current i','LineWidth',2,'fontsize',14)
%title('full state observer response i_0 = 0.75 amper')
%legend('plant','observer')
%figure(6)
%plot(ans.tout,ans.thd_p,'b','LineWidth',2)
%hold on
%plot(ans.tout,ans.thd_ob,'r','LineWidth',2)
%xlabel('Time','LineWidth',2,'fontsize',14)
%ylabel('rate of \theta (radian/second)','LineWidth',2,'fontsize',14)
%title('full state observer response rate of \theta_0 = 0.5 rad/s')
%legend('plant','observer')
%figure(7)
%plot(ans.tout,ans.phid_p,'b','LineWidth',2)
%hold on
%plot(ans.tout,ans.phid_ob,'r','LineWidth',2)
%xlabel('Time','LineWidth',2,'fontsize',14)
%ylabel('rate of \phi (radian/second)','LineWidth',2,'fontsize',14)
%title('full state observer response rate of \phi_0 = 0.5 rad/s')
%legend('plant','observer') 
 