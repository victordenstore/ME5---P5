% The purpose of this model is to get a tool to analyze an immersion type transudcer which is used to manipulate the position of air bubbles in oil. 
% The model is based on a model created by Sittig and utilizes derivations of it in the form of transfer functions. An expression for the pressure field
% was derived from the 1D wave equation. Eventually the force on the particles was calculated by utilization of Gor'kov's potential.
% 
% The authors of this script are Victor Dahl, Phillip Kjaer, Wenzel Neumann
% and Jeppe Ottense


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
% dfr=0.0000042;                       % Random variable to establish acoustic backing material relationship
% Z_b=dfr*w;                           % Acoustic impedance of the backing plate (this is a function of frequency)  
k=w/v_o;                           % wave number for the peizoelectric plate - found from NDE book


%% Acoustic impedance calculations
rho_4=7860;                          % Density of plain carbon steel 
c_3=3230;                            % Speed of wave in plain carbon steel
Z_b=rho_4*c_3*S;                     % Acoustic impedance of the backing plate

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


%% Matricies used in the Sittig Model

TAe_matrix=[1/n n/(1i*w*C_o); -1i*w*C_o 0];                                                                              %Transducer electrical matrix

TAa_matrix=(1/(Z_b-1i*Z_o*tan(k*d/2)))*[Z_b+1i*Z_o*cot(k*d) (Z_o)^2+1i*Z_o*Z_b*cot(k*d); 1 Z_b-2*1i*Z_o*tan(k*d/2)];     %Transducer acoustic matrix

TAl_matrix=[cos(k_2*d_2) -1i*Z_m*sin(k_2*d_2); (-1i*sin(k_2*d_2))/Z_m cos(k_2*d_2)];                                     %Transducer coating material matrix



TA_matrix=TAe_matrix*TAa_matrix*TAl_matrix;                                                                              %Sittig model matrix

%% equations for Sensitvities

S_vI=1/(Z_r*TA_matrix(2,1)+TA_matrix(2,2));                                   %Sensitivity for vI.

Z_in=(Z_r*TA_matrix(1,1)+TA_matrix(1,2))/(Z_r*TA_matrix(2,1)+TA_matrix(2,2)); %Electrical impedance of the transducer

S_FV=Z_r*S_vI/Z_in;                                                           % Sensitivity FV

%% Pressure at the transducer surface

V_M = 150;                                                                    % Input voltage
lambda =c/f;                                                                  % wavelength of the pressure field
F_surface = V_M * S_FV;                                                       % Total force on trasnducer surface
P_surface = F_surface/(S_a);                                                  % pressure on the transducer surface

%% Translation to the 1D pressure field

damp_coeff =(2*Kin_Vis*w^2)/(3*c^3);                                          % Damping coefficient for attenuation within the liquid
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

t = linspace(0,10*pi/w,1000);                                               % time duration during which the pressure field
                                                                            % is monitored.
x = linspace(0,5*0.0013,1000);                                              % Distance over which the pressure field is monitored.

% Here the P_field matrix is filled with the pressure amplitudes and phase
% shifts at a certain location in space and moment in time.
 for n = 1:1000
     for m = 1:1000
     P_field(n,m) = P_surface*(R1+1) * exp(1i*w*t(m)-1i*x(n)*2*pi/lambda) ...
         - R2*P_surface * exp(1i*w*t(m)+1i*(x(n))*2*pi/lambda);
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


 
%% Constants used to calculate the acoustic radiation force

rho_p = 1.225;                                                                  %Air density [kg/m^3]
r = 80*10^-6;                                                                   %Bubble radius [m]
Vp = (4*pi/3)*r^3;                                                              %Bubble volume [m^3]
Mp = rho_p * Vp;                                                                %Bubble mass [kg]

% Oil of type ISO VG32 is used as the fluid
rho_l = 868;                                                                    % Oil density [kg/m^3]
mu = 32*10^-6*rho_l;                                                            % Oil Viscosity [kg/m*s]

% Besides the acoustic radiation force, acoording to Stoke's Law, there is
% also a drag force acting on the particle, which is given by f_D = B*v

B = 6*pi*r*mu;                                                                  % coefficient (gain) of Stoke's Law

% The acoustic radiation force is dependent on the compressibility factor
% f_1 and the density factor f_2

K_p = 1/(rho_p*330^2);                                                      % Compressibility of the air bubble
K_l = 1/(rho_l*c^2);                                                        % Compressibility of the oil
f_1=1-(K_p/K_l);                                                            % Compressibility factor
f_2=2*(rho_p-rho_l)/(2*rho_p+rho_l);                                        % Density factor

% The acoustic contrast factor determines whether the particles are forced
% towards the pressure nodes or antinodes. A negative contrast factor means
% that particles are forced towards the antinodes, whils a positive
% indicates movement towards the nodes.

acf = (5*rho_p-2*rho_l)/(2*rho_p+rho_l) - (K_p/K_l);                        % The acoustic contrast factor

% To calculate the acoustic radiation force, the time averaged pressure is
% required. For a standing wave it suffices to take only one period into
% account.

T = 1/f;                                                                    % One period of the pressure field

%% Plotting the results obtained from the acoustic radiation force

velocities = zeros(10,1);                                                   % A vector is initialized to store the values of  
                                                                            % velocities at different locations in x direction.

for i = 1:length(velocities)
    x = linspace(0.0013,0.0026,length(velocities));                         % The distance over which the velocity is monitored 
    x = x(i);                                                                        
v = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B);    % Call upon the function to obtain the velocity.
vnew = subs(v);                                                             % Constants are substituted for the symbolic variables.
velocities(i) = double(vnew);                                               % the velocities are stored as scalars.
end

x = linspace(0.0013,0.0026,length(velocities));
% Here the resultant velocity is plotted for different locations x away
% from the transducer surface.
figure
plot(x',velocities)
xlabel('distance [m]'); ylabel('velocity [m/s]'); title('Particle velocities');


%% Here I tried to plot the different contributions to the acoustic force similar to figure 9 of the Trujillo paper

x = linspace(0.0013,0.0026,length(velocities));                             % Monitored distance, set to one wavelength, can be changed. Called again x
                                                                            % because it is the same variable in the Gorkov function but the domain can be 
                                                                            % changed.
% Here the outputs of the function are retrieved and converted from 
% symbols to numbers.
[v,p_mean_square,v_in,v_mean_square,Uac,Fac]  = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B); 
p_ms = double(subs(p_mean_square));
v_ms = double(subs(v_mean_square));
Uac = double(subs(Uac));
Fac = double(subs(Fac));

p_ms_vc = [p_ms; p_ms; p_ms; p_ms; p_ms; p_ms; p_ms; p_ms; p_ms; p_ms];     % Matrices for plotting. Constructed 
v_ms_vc = [v_ms; v_ms; v_ms; v_ms; v_ms; v_ms; v_ms; v_ms; v_ms; v_ms];     % in the same way as for the pressure
Uac_vc = [Uac; Uac; Uac; Uac; Uac; Uac; Uac; Uac; Uac; Uac];                % field plot, only this time there is 
Fac_vc = [Fac; Fac; Fac; Fac; Fac; Fac; Fac; Fac; Fac; Fac];                % no change in values through time, as 
                                                                            % they are time averages.
time_averaged_vec = linspace(0,T,length(velocities));                       % Since all are time-averaged, I just made a random duration 
                                                                            % of one period,
                                                                            % which I use to be
                                                                            % able to make a
                                                                            % contour plot.
% Plots according to figure 9 of Trujillo                                                                            
                                                                            
figure                                                      
tiledlayout(4,1); nexttile;
contourf(x,time_averaged_vec,real(p_ms_vc), 14); colormap jet; colorbar;
xlabel('distance [m]'); ylabel('time [s]'); title('<p^2>');

nexttile;
contourf(x,time_averaged_vec,real(v_ms_vc), 14); colormap jet; colorbar;
xlabel('distance [m]'); ylabel('time [s]'); title('<v^2>');

nexttile;
contourf(x,time_averaged_vec,real(Uac_vc), 14); colormap jet; colorbar;
xlabel('distance [m]'); ylabel('time [s]'); title('Uac');

nexttile;
contourf(x,time_averaged_vec,real(Fac_vc), 14); colormap jet; colorbar;
xlabel('distance [m]'); ylabel('time [s]'); title('Fac');








