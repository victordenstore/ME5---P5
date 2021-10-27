 clear all; close all; clc;

%% Initial Variables
f = 1*10^6;                      % system frequency in Hz
w=2*pi*f;                        % System frequency in rad/s

%% Piezoelectric constants

r_p=2.5*10^-3;                     % Radius of the piezoelectric plate - Found from CAD model
rho=7.8*10^3;                      % Piezoelectric plate density - Found from data sheet
d=2*10^-3;                         % Piezoelectric plate thickness - found from data sheet
S=pi*r_p^2;                        % Piezoelectric surface area - area of a circle due to circular peizoelectric disk face

c_33=16.6*10^10;                   % Elastic constant of the plate - found from source (victor found - please put it in when editing)
d_33=265*10^-12;                   % Piezoelectric charge coefficient - found from source -||-
s_33=14.2*10^-12;                  % Elastic compliance coefficient - found from source -||-
epsilon_0=8.854*10^-12;            % Vacuum permitivity - found from source -||-
epsilon_33=1200*epsilon_0;         % relativity Permativity - found from source -||-
Beta_33= epsilon_0/epsilon_33;     % the dielectric impermeability of the plate at constant strain - found from source -||-
h_33= d_33/(s_33*epsilon_33);      % Piezoelectric stiffness constant for the plate - found from source -||-
v_o=sqrt(c_33/rho);                % compressional wave speed from piezoelectric - found from source -||-

%%TAe calculated inputs
C_o=S/(Beta_33*d);                 % the clamped capacitance of the plate - Found from NDE book
n=h_33*C_o;                        % A given constant

%TAa calculated inputs
x=0.0000042;                       % Random variable to establish acoustic backing material relationship
Z_b=x*w;                           % Acoustic impedance of the backing plate (this is a function of frequency)  
k=w/v_o;                           % wave number for the peizoelectric plate - found from NDE book


%% Acoustic impedance calculations
Z_o=rho*v_o*S;                       % plane wave acoustic impedance of the piezoelectric plate

S_a=S;                               % effective face area of the transducer
rho_2=857;                           % density of the fluid ISO VG 32
c=1300;                              % compresisonal wave speed in fluid 
Z_r=S_a*c*rho_2;                     % Acoustic radiation impedance 

c_2=3230;                            % Speed of sound in steel
rho_3=7850;                          % Density of steel
k_2=w/c_2;                           % wave number of the coating material of the transducer source
d_2=0.5*10^-3;                       % thickness of the coating material of the transducer source
S_m=(4.5*10^-3)^2*pi;                % Face area of the coating material fo the transducer soruce
Z_m=rho_3*c_2*S_m;                   % The acoustic impedance of the coating material of the transducer source


Kin_Vis=32*10^-6;                % Kinematic Viscosity of the fluid (ISO VG 32)
Vis=Kin_Vis*rho_2;               % Dynamic viscosity of the fluid
R=((c*rho_2)-(c_2*rho_3))/((c*rho_2)+(c_2*rho_3)); %Reflection coefficient 
L=0.04485;                                         %Lentgh of chamber


%% Matricies
TAe_matrix=[1/n n/(1i*w*C_o); -1i*w*C_o 0];                                                                              %Transducer electrical matrix

TAa_matrix=(1/(Z_b-1i*Z_o*tan(k*d/2)))*[Z_b+1i*Z_o*cot(k*d) (Z_o)^2+1i*Z_o*Z_b*cot(k*d); 1 Z_b-2*1i*Z_o*tan(k*d/2)];     %Transducer acoustic matrix

TAl_matrix=[cos(k_2*d_2) -1i*Z_m*sin(k_2*d_2); (-1i*sin(k_2*d_2))/Z_m cos(k_2*d_2)];                                     %Transducer coating material matrix



TA_matrix=TAe_matrix*TAa_matrix*TAl_matrix;                                                                              %Sittig model matrix

%% equations for Sensitvities

S_vI=1/(Z_r*TA_matrix(2,1)+TA_matrix(2,2)); %Sensitivity for vI.

Z_in=(Z_r*TA_matrix(1,1)+TA_matrix(1,2))/(Z_r*TA_matrix(2,1)+TA_matrix(2,2)); %Electrical impedance of the transducer

S_FV=Z_r*S_vI/Z_in; %Sensitivity FV

%% Graphing

t = linspace(0,10*pi/w,1000);
x = linspace(0,L,1000);
alfa = -40;
V_in = zeros(1000,1);
%F = zeros(1000,1);
%P = zeros(1000,1);
%expression = zeros(1000,1000);

%% Solving for Pressure at transducer face
S_FV = abs(S_FV)*exp(1i*atan((imag(S_FV)/real(S_FV))))
V_M = 150;
lambda =c/f;
F_surface = V_M * S_FV;
P_surface = F_surface/S_a;

%% Solving for Pressure field

damp_coeff =(2*Kin_Vis*w^2)/(3*c^3);
exp1 = zeros(1000,1000);
shit = zeros(1000,1);

 for to = 1:1000
     for gi = 1:1000
     exp1(to,gi) = P_surface*(4.2736+1) * exp(1i*w*t(gi)-1i*x(to)*2*pi/lambda)*exp(- damp_coeff * x(to))-4.8295*P_surface * exp(1i*w*t(gi)-1i*(-x(to))*2*pi/lambda)*exp(- damp_coeff * (abs(x(to)-L)));
     
     end
 end




contourf(t, x, real(exp1), 14); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); 

figure
time_vecs = linspace(0,10*pi/w,1000);

% for phil = 1:1000
meshc(t, x, real(exp1)); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); zlabel('pressure [Pa]')

 %hold on
 %end

%% Constants for Simulink matlab code:
rho_p = 1.225; %Air density kg/m^3
r = 1000; %Bubble radius m
Vp = 4*pi/3*r^3; %Bubble volume m^3
Mp = rho_p * Vp; %Bubble mass kg

rho_l = 868; %Oil density kg/m^3
mu = 32*10^-6*rho_l; %ISO VG32 Viscosity kg/m*s
B = 6*pi*r*mu;

K_p = 1; %Compressibility factor of the particle
K_l = 1/(1.8*10^4); %Of the liquid

f_1=1-(K_p/K_l); %f_1 of the Gorkov eq.
f_2=2*(rho_p-rho_l)/(2*rho_p+rho_l); %f_2 of the Gorkov eq.
P = P_surface;


% a = zeros(10,1);
% x = zeros(11,1);
% t = linspace(0,11*2*pi/w,11)';
v = zeros(11,1);
T = 1/f;
x = 0; %5cm, intial bubble placement
v = 25;       % intital bubble velocity
t = 0;       % initial time
% a = EOM_Particle(Mp,Gorkov(c, rho_l, r, f_1, f_2,w,T,x(1),damp_coeff,P),B,v(1));
acoustic_contrast_factor = (5*rho_p-2*rho_l)/(2*rho_p+rho_l) - (K_p/K_l)
Eac = P_surface*2 / (4*rho_l*c^2);


%a = EOM_Particle(Mp,Gorkov(c, rho_l, r, f_1, f_2,w,T,x,damp_coeff,P_surface),B,v)
a = EOM_Particle(Mp,Gorkov(c, rho_l, r, f_1, f_2,w,x,P_surface,k),B,v)
t1 = 0;
t2 = T;
v= integral(@(t) (a),t1,t2,'ArrayValued',true)
v = abs(v)



% for i = 1:10
% a(i) = EOM_Particle(Mp,Gorkov(c, rho_l, r, f_1, f_2,w,T,x(i),damp_coeff,P),B,v(i));
% a_new = a(i);
% t1 = t(i);
% t2 = t(i+1);
% v(i+1) = integral(@(t) (a_new),t1,t2,'ArrayValued',true);
% v_new = v(i+1);
% x(i+1) = integral(@(t) (v_new),t1,t2,'ArrayValued',true);
% end

% figure
% plot(t,x)
% 



% %% Open Simulink Model 
% open('SimulinkModel');
% 
% 
% %% Simulation
% Simulation = sim('SimulinkModel');

