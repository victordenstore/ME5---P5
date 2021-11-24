close all;
clear all;
clc;
%% Parameter initialization

j = 1i;                                                                     % imaginary value
f = 4.30*10^5;                                                               % operating frequency [Hz]
omega = 2*pi*f;                                                             % operating frequency [rad/s] 
c_33 = 16.6*10^10;                                                          % Elastic compliance coefficient
c_n = c_33;                                                                 % Comes from the multilayer paper.
rho_n = 7.8*10^3;                                                           % density of the piezoelectric disc.
v_n = sqrt(c_n/rho_n);                                                      % wave speed of compressional waves in the piezoelectric plate
c_mn = v_n;                                                                 % sound velocity.
d_n = 2*10^-3;                                                              % thickness of the piezoelectric plate.
l_mn = d_n;                                                                 % thickness according to Sittig schwanzel.
f_mn = c_mn/(2*l_mn);                                                       % some specific frequency I guess.
%f_backing = c_backing/(2*l_backing);


gamma = omega * d_n/v_n;                                                    % 
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
s_33 = 14.2*10^-12;
d_33 = 265*10^-12;
h_33 = d_33/(s_33*eps_33);
k_n_squared = h_33^2 *eps_33/c_33;
k_squared = k_n_squared;
s = k_squared * sin(gamma)/gamma;
%gamma_mn = [pi*f/f_mn pi*f/f_backing];
r = 2.5*10^-3; 
S = pi*r^2;
C_0= S*eps_33/(d_n);
Z_0 = rho_n*v_n*S;
omega_mn = 2*pi*f_mn;
Z_mn = Z_0;
phi = sqrt(k_squared)*sqrt(omega_mn*c_mn*Z_mn/pi);
c = k_squared * (1-cos(gamma))/gamma;
A_mn = S;
Q_mn = rho_n;
Z_mn = A_mn * Q_mn; 

rho_Back = 7850;
c_Back = 3230;
r_Back = 2.5*10^-3;
S_Back = pi*r_Back^2;
Z_Back = rho_Back*c_Back*S_Back;  
Z_NB = Z_Back;       % not sure about this! As NB is the final layer, but front or back??

c_fluid = 1300;
Kin_Vis=32*10^-6; 




%% 4 port model matrix - Sittig


T_11 = (cos(gamma) - s)/(1-s);
T_12 = (j*Z_0*(sin(gamma) - 2*c))/(1-s);
T_13 = -((cos(gamma) - 1)*phi)/(1-s);
T_14 = 0;
T_21 = (j*sin(gamma))/(Z_0*(1-s));
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

Tn=[T_11 T_12 T_13 T_14; ...
    T_21 T_22 T_23 T_24; ...
    T_31 T_32 T_33 T_34; ...
    T_41 T_42 T_43 T_44];
T=(Tn)^2;
T_11=Tn(1,1);
T_12=Tn(1,2);
T_13=Tn(1,3);
T_14=Tn(1,4);
T_21=Tn(2,1);
T_22=Tn(2,2);
T_23=Tn(2,3);
T_24=Tn(2,4);
T_31=Tn(3,1);
T_32=Tn(3,2);
T_33=Tn(3,3);
T_34=Tn(3,4);
T_41=Tn(4,1);
T_42=Tn(4,2);
T_43=Tn(4,3);
T_44=Tn(4,4);
%% Intermediate layers (ex. glue layers, copper electrode, etceteraness)

r_copper = 2.5*10^-3;
A_copper = pi*r_copper^2;
rho_copper = 8920;
c_copper = 3570;
d_copper = 0.5*10^-3;
v_copper = sqrt(c_copper/rho_copper)

Z_copper = A_copper*rho_copper*c_copper;
gamma_copper = omega * d_copper/v_copper;

A_copper = [cos(gamma_copper) j*Z_copper*sin(gamma_copper); ...
    (j*sin(gamma_copper))/Z_copper cos(gamma_copper)];

rho_FP = 7850;
r_FP = 4.5*10^-3;
c_FP = 3230;
S_FP = pi*r_FP^2;
d_FP = 0.5*10^-3;
v_FP = sqrt(c_FP/rho_FP)

Z_FP = rho_FP*c_FP*S_FP;
gamma_FP = omega*d_FP/v_FP;

A_FP = [cos(gamma_FP) j*Z_FP*sin(gamma_FP); ...
    (j*sin(gamma_FP))/Z_FP cos(gamma_FP)];

rho_H = 7860;
r_H = 2.5*10^-3;
c_H = 3230;
S_H = pi*r_H^2;
d_H = 3.25*10^-3;
v_H = sqrt(c_H/rho_H);
gamma_H = omega*d_H/v_H;
Z_H = rho_H*c_H*S_H;
A_H = [cos(gamma_H) j*Z_H*sin(gamma_H); ...
    (j*sin(gamma_H))/Z_H cos(gamma_H)];

A_m = A_copper*A_FP*A_H;
A_b = A_m(1,1);
B_b = A_m(1,2);
C_b = A_m(2,1);
D_b = A_m(2,2);

Z_b = (A_b*Z_NB + B_b)/(C_b*Z_NB + D_b);

%% 2x2 matrices for parallel connection of piezoelectric layers

A_c = T_31 - T_33*(T_21*Z_b + T_11)/(T_23*Z_b + T_13);
B_c = T_32 - T_33*(T_22*Z_b + T_12)/(T_23*Z_b + T_13);
C_c = T_41 - T_43*(T_21*Z_b + T_11)/(T_23*Z_b + T_13);
D_c = T_42 - T_43*(T_22*Z_b + T_12)/(T_23*Z_b + T_13);

A_C = [A_c B_c; C_c D_c];

%% Final 2x2 matrix of Piezo stack 
A_ps = A_C*A_m;
A = A_ps(1,1);
B = A_ps(1,2);
C = A_ps(2,1);
D = A_ps(2,2);


%% Sensitivity Functions

S_a=S;                               % effective face area of the transducer
rho_2=857;                           % density of the fluid ISO VG 32
c=1300;                              % compresisonal wave speed in fluid 
Z_r=S_a*c*rho_2;                     % Acoustic radiation impedance 


Z_E = (A*Z_H + B)/(C*Z_H + D);            % electrical input impedance
V_IL = Z_H/(A*Z_H + B);                   % voltage transfer ratio between input voltage and output force



V_M = 150;                                                                    % Input voltage
lambda =c_fluid/f;                                                                  % wavelength of the pressure field
F_surface = V_M * V_IL;                                                       % Total force on trasnducer surface
P_surface = F_surface/S;                                                  % pressure on the transducer surface

%% Translation to the 1D pressure field

damp_coeff =(2*Kin_Vis*omega^2)/(3*c_fluid^3);                                          % Damping coefficient for attenuation within the liquid
P_field = zeros(1000,1000);                                                   % A matrix is initialized for storing the (complex)
                                                                              % pressure values varying with time and distance.     
% A standing wave is produced by attenuation of the travelling wave and its 
% reflections resulting from rebounding at the boundaries (walls). During
% the reflection at the boundaries, energy is lost which can be represented
% by a reflection coefficient, which takes the number of reflections until
% total attenuation into account. There is distinguished between the
% reflection coefficient due to impact with the right and left wall.

R1 = 4.2736;                                                                % Left wall reflection coefficient.
R2 = 4.8295;                                                                % Right wall reflection coefficient.

t = linspace(0,10*pi/omega,1000);                                               % time duration during which the pressure field
                                                                            % is monitored.
x = linspace(0,5*0.0013,1000);                                              % Distance over which the pressure field is monitored.

% Here the P_field matrix is filled with the pressure amplitudes and phase
% shifts at a certain location in space and moment in time.
 for n = 1:1000
     for m = 1:1000
     P_field(n,m) = P_surface*(R1+1) * exp(1i*omega*t(m)-1i*x(n)*2*pi/lambda) ...
         - R2*P_surface * exp(1i*omega*t(m)+1i*(x(n))*2*pi/lambda);
      end
 end

%% Contourplots of the standing pressure wave

% 2D plot which shows the pressure minima and maxima of the pressure field
figure
contourf(t,x, real(P_field), 14); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); title('Contourplot of the 1D pressurefield');

%% Surface plots of the standing pressure wave

% 3D animation of the pressure field
figure
meshc(t,x, real(P_field)); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); zlabel('pressure [Pa]'); title('Surface plot of the 1D pressure field');









    
    