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
J_H=0.001; % Housing's inertia
J_m=0.032; % motor's shaft inertia
j=9.7*10^-7; % gearbox's output shaft inertia
R_m=1.6; % DC motor's Electric resistance
L_m=0.145; % DC motor's inductance
K_m=0.0109; % motor's constant
N=14; % Gearbox transfer ratio
g=9.8; % gravity acceleration
N=14; % Gearbox transfer ratio
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




%%simulation
load_system('linear_model_advanced_control_project')
        sim('linear_model_advanced_control_project')

        
        
 figure(1)
 plot(ans.tout,ans.tl,'k')

 xlabel('Time','LineWidth',2,'fontsize',14)
 ylabel('\theta','LineWidth',2,'fontsize',14)
 
 legend
