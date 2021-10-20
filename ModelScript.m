 clear all; close all; clc;
x=0.0000042;

%% Initial Variables
w=2*pi*10^6;                     % System frequency in rad/s
a=2.5*10^-3;                     % Radius of the piezoelectric plate
rho=7.8*10^3;                    % Piezoelectric plate density
c_33=16.6*10^10;                 % Elastic constant of the plate
d=2*10^-3;                       % plate thickness
Z_b=x*w;                         % Acoustic impedance of the backing plate (this is a function of frequency)
d_33=265*10^-12;                 % Piezoelectric charge coefficient
s_33=14.2*10^-12;                % Elastic compliance coefficient
epsilon_0=8.854*10^-12;          % Vacuum permitivity
epsilon_33=1200*epsilon_0;       % relativity Permativity
Beta_33= epsilon_0/epsilon_33;   % the dielectric impermeability of the plate at constant strain,
c=1300;                          % compresisonal wave speed in fluid
f=w/(2*pi);                      % System Angular frequency 1/s
h_33= d_33/(s_33*epsilon_33);    % Piezoelectric stiffness constant for the plate
v_o=sqrt(c_33/rho);              % compressional wave speed from piezoelectric
k=w/v_o;                         % wave number for the peizoelectric plate
S=pi*a^2;                        % Piezoelectric surface area
C_o=S/(Beta_33*d);               % the clamped capacitance of the plate
n=h_33*C_o;                      % A given constant
S_a=S;                           % effective face area of the transducer
rho_2=857;                       % density of the fluid
Kin_Vis=32*10^-6;                % Kinematic Viscosity of the fluid (ISO VG 32)
Vis=Kin_Vis*rho_2;               % Dynamic viscosity of the fluid

Z_o=rho*v_o*S;                   % plane wave acoustic impedance of the piezoelectric plate
c_2=3230;                        %Speed of sound in steel
rho_3=7850;                      %Density of steel
R=((c*rho_2)-(c_2*rho_3))/((c*rho_2)+(c_2*rho_3)); %Reflection coefficient 
L=0.00448;                        %Lentgh of chamber

%% Matricies
TAe_matrix=[1/n n/(1i*w*C_o); -1i*w*C_o 0];          %Transducer electrical matrix
TAa_matrix=(1/(Z_b-1i*Z_o*tan(k*d/2)))*[Z_b+1i*Z_o*cot(k*d) (Z_o)^2+1i*Z_o*Z_b*cot(k*d); 1 Z_b-2*1i*Z_o*tan(k*d/2)];     %Transducer acoustic matrix

TA_matrix=TAe_matrix*TAa_matrix; %Sittig model matrix

%% equations for acoustic Impedance
Z_r=S_a*c*rho_2;         % Acoustic radiation impedance - just in case we want to put in the Ka condition later.

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

V_M = 150;
lambda =c/f;
P = S_FV * V_M/S;

%% Solving for Pressure field

damp_coeff =(2*Kin_Vis*w^2)/(3*c^3);
exp1 = zeros(1000,1000);

 for to = 1:1000
     for gi = 1:1000
     exp1(to,gi) = P * exp(1i*w*t(gi)-1i*x(to)*2*pi/lambda - damp_coeff * x(to))+exp1(to,gi)+R*P * exp(1i*w*t(gi)+1i*pi-1i*(-x(to))*2*pi/lambda - damp_coeff * (x(to)+L));
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
% rho_p = 1.225; %Air density kg/m^3
% r = 1*10^-3; %Bubble radius m
% Vp = 4*pi/3*r^3; %Bubble volume m^3
% Mp = rho_p * Vp; %Bubble mass kg
% 
% rho_l = 868; %Oil density kg/m^3
% mu = 32*10^-6*rho_l; %ISO VG32 Viscosity kg/m*s
% B = 6*pi*r*mu;
% 
% K_p = 1; %Compressibility factor of the particle
% K_l = 1/(1.8*10^4); %Of the liquid
% 
% f_1=1-(K_p/K_l); %f_1 of the Gorkov eq.
% f_2=2*(rho_p-rho_l)/(2*rho_p+rho_l); %f_2 of the Gorkov eq.


%a = zeros(10,1);
%x = zeros(11,1);
%t = linspace(0,11*2*pi/w,11)';
%v = zeros(11,1);
%x(1) = 5*10^-2; %5cm, intial bubble placement
%v(1) = 0;       % intital bubble velocity
%t(1) = 0;       % initial time

%for i = 1:10
%a(i) = EOM_Particle(Mp,Gorkov(c, rho_l, r, f_1, f_2,w,t(i),x(i),damp_coeff,P,lambda),B,v(i));
%a_new = a(i);
%t1 = t(i);
%t2 = t(i+1);
%v(i+1) = integral(@(t) (a_new),t1,t2,'ArrayValued',true);
%v_new = v(i+1);
%x(i+1) = integral(@(t) (v_new),t1,t2,'ArrayValued',true);
%end


%plot(t,x)




% %% Open Simulink Model 
% open('SimulinkModel');
% 
% 
% %% Simulation
% Simulation = sim('SimulinkModel');

