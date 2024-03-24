clc;clear
%%
A_ee = [0 0 0 0; 0 0 1 0; 5.5413 0 0 -3.9; 0 0 -0.2255 -11.03];
syms l1 l2 l3 l4  i u x_e1 x_e2 x_e3 x_e4 s
C_s = 1;
A_se = [-16.26 0 0 0.25];
A_es = [-16.26 0 5.5 0]';
L = [l1 l2 l3 l4]';
A_ss = 0;
B_e = [0 0 0 0]';
B_s = 6.8966;
x_e = [x_e1 x_e2 x_e3 x_e4]';
cl = (A_ee - L*C_s*A_se)*x_e + (A_es - L*A_ss)*i + (B_e - L*B_s)*u;
cll = det(s*eye(4,4)-A_ee+L*A_se)
r=rank(obsv(A_ee,A_se))

%% gain designing
tt = 5*[-8 -1.23 -7 -3];
 l = place(A_ee,A_se',tt)


%l_1 = -29/0.2255;
%l_2 = -233.96;
%l_4 = -307;
%l_3 = 100;