close all;
clear all;
clc;


%% Parameter initialization

j = 1i;
f = 4.3*10^5;
omega = 2*pi*f; 

v_n = sqrt(c_n/rho_n);
c_mn = v_n;
d_n = 2*10^-3;
l_mn = d_n;
f_mn = c_mn/(2*l_mn);
c_backing = 
f_backing = c_backing/(2*l_backing);
rho_n = 7.8*10^3; 
c_33 = 16.6*10^10;
c_n = c_33;


gamma = omega * d_n/v_n;
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
s_33 = 14.2*10^-12;
d_33 = 265*10^-12;
h_33 = d_33/(s_33*eps_33);
k = sqrt(h_33^2 *eps_33/c_33);
s = k^2 * sin(gamma)/gamma;
gamma_mn = [pi*f/f_mn pi*f/f_backing];
r = 2.5*10^-3; 
S = pi*r^2;
Z_0 = rho_n*v_n*S;
omega_mn = 2*pi*f_mn;
Z_mn = Z_0;
phi = k*sqrt(omega_mn*c_mn*Z_mn/pi);
c = k^2 * (1-cos(gamma))/gamma;
A_mn = S;
Q_mn = rho_n;
Z_mn = A_mn * Q_mn; 
rho_NB = 7850;
c_NB = 3230;
S_NB = pi*(4.5*10^-3)^2
Z_NB = rho_NB*c_NB*S_NB;  % not sure about this one!

%% 4 port model matrix


T_11 = (cos(gamma) - s)/(1-s);
T_12 = (j*Z_0*(sin(gamma) - 2*c))/(1-s);
T_13 = -((cos(gamma) - 2*c)*phi)/(1-s);
T_14 = 0;
T_21 = (j*sin(gamma))/(Z_0(1-s));
T_22 = (cos(gamma) - s)/(1-s);
T_23 = -(j*phi*sin(gamma))/(1-s);
T_24 = 0;
T_31 = 0;
T_32 = 0;
T_33 = 1;
T_34 = 0;
T_41 = -(j*sin(gamma))/(Z_0*(1-s)) * phi;
T_42 = -((cos(gamma) - 1)*phi)/(1-s);
T_43 = (j*omega*C_0)/(1-s);
T_44 = 1;

%% Intermediate layers

A_mn = [cos(gamma_mn) j*Z_mn*sin(gamma_mn); ...
    (j*sin(gamma_mn))/Z_mn cos(gamma_mn)];


Z_b = (A_b*Z_NB + B_b)/(C_b*Z_NB + D_b);

%% 2x2 matrices for parallel connection of piezoelectric layers

A_c = T_31 - T_33*(T_21*Z_b + T_11)/(T_23*Z_b + T_13);
B_c = T_32 - T_33*(T_22*Z_b + T_12)/(T_23*Z_b + T_13);
C_c = T_41 - T_43*(T_21*Z_b + T_11)/(T_23*Z_b + T_13);
D_c = T_42 - T_43*(T_22*Z_b + T_12)/(T_23*Z_b + T_13);



    
    