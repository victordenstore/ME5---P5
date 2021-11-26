% This script is an intent to accurately model a piezoelectric transuducer,
% consisting of multiple intermediate piezoelectric inactive layers. The
% modelling is done following the approach as described by Sittig in his
% paper: 'Effects of Bonding and Electrode Layers on the  Transmission     
% Parameters of Piezoelectric Transducers Used in Ultrasonic Digital Delay 
% Lines'. At first a model will be made consisting of one active
% piezoelectric plate, a backing side consisting of backing material and a
% copper electrode, and the transmission side consisting of a copper
% electrode, a front plate and housing material. The models transfer
% function from voltage input to force output, the impedance frequency  
% response function and the admittance frequency response function will be
% analyzed and if satisfying results are obtained, there will be extended
% upon the model by inclusion of a second piezoelectric active layer in
% parallel with the first layer and by taking glue layers into
% consideration.

%% close all former work
clear all; close all; clc;

%% Parameter initialization

f = 0.43*10^6;      
omg = 2*pi*f;
j = 1i;
%% Setting up the backing material matrix

% First the layer closest to the piezoelectric layer, hence the copper
% electrode.

% the material parameters

rho_Cu = 8960;  % Pure! copper density.
c_Cu = 3570;    % pure copper speed of sound.
l_Cu = 0.5*10^-3;
r_Cu = 2.5*10^-3;
A_Cu = pi*r_Cu^2;
Z_0_Cu = rho_Cu * c_Cu * A_Cu;
f_0_Cu = c_Cu/(2*l_Cu);
gamma_Cu = pi*f/f_0_Cu;

% The matrix

Back_mat_Cu = [cos(gamma_Cu) j*Z_0_Cu*sin(gamma_Cu); ...
               (j*sin(gamma_Cu))/Z_0_Cu cos(gamma_Cu)];
           
           
% Similarily the matrix is set up for the thick, cone shaped backing
% material

% the material parameters

rho_FeC = 8000;  % Pure! copper density.
c_FeC = 5800;    % pure copper speed of sound.
r_FeC = 2.5*10^-3;
A_FeC = pi*r_FeC^2;
Z_0_FeC = rho_FeC * c_FeC * A_FeC;
f_0_FeC = c_FeC/(2);                     % As the thickness of the backing is assumed to approach infinity, it is omitted. Not sure about this! Check!!
gamma_FeC = pi*f/f_0_FeC;

% the matrix

Back_mat_FeC = [cos(gamma_FeC) j*Z_0_FeC*sin(gamma_FeC); ...
                (j*sin(gamma_FeC))/Z_0_FeC cos(gamma_FeC)];
            
% Setting up the matrices for the additional glue layers between the back face of the piezoelectric element and the 
% copper electrode, and between the copper electrode and the backing plate.

rho_glue = 960;  
c_glue = 3570;    % unsure about this!!
l_glue = 0.01*10^-3;    % to be varied.
r_glue = 2.5*10^-3;
A_glue = pi*r_Cu^2;
Z_0_glue = rho_Cu * c_Cu * A_Cu;
f_0_glue = c_Cu/(2*l_Cu);
gamma_glue = pi*f/f_0_Cu;

glue_mat = [cos(gamma_glue) j*Z_0_glue*sin(gamma_glue); ...
               (j*sin(gamma_glue))/Z_0_glue cos(gamma_glue)];
            
% The matrices are multiplied and the total impedance Z_b of the backing
% layers is obtained.

Back_mat = glue_mat*Back_mat_Cu*glue_mat*Back_mat_FeC;
A_b = Back_mat(1,1);
B_b = Back_mat(1,2);
C_b = Back_mat(2,1);
D_b = Back_mat(2,2);

Z_0b = Z_0_FeC;
Z_b = (A_b*Z_0b + B_b)/(C_b*Z_0b + D_b);    % Here Z_0b is the impedance of the final backing layer only!

%% Setting up the transmission layer matrix

% The intermediate layers in the transmission direction are obtained in the
% same way. Again starting with the intermediate layer closest to the
% peizoelectric active layer.

% The material parameters were already given before!

Trans_mat_Cu = [cos(gamma_Cu) j*Z_0_Cu*sin(gamma_Cu); ...
               (j*sin(gamma_Cu))/Z_0_Cu cos(gamma_Cu)];
% For now the glue layer is omitted and the matrix is setup for the front
% plate only. Which is made out of Steel, just as the backing material.

% material properties of steel front plate.
r_FP = 4.5*10^-3;
A_FP = pi*r_FP^2;
Z_0_FP = rho_FeC * c_FeC * A_FP;
l_FeC = 0.5*10^-3;
f_0_FP = c_FeC/(2*l_FeC);                    
gamma_FeC = pi*f/f_0_FP;

Trans_mat_FP = [cos(gamma_FeC) j*Z_0_FP*sin(gamma_FeC); ...
                (j*sin(gamma_FeC))/Z_0_FP cos(gamma_FeC)];

% Multiplication of the transmission piezoelectric inactive layers gives:

Trans_mat = glue_mat*Trans_mat_Cu*glue_mat*Trans_mat_FP;
                   
%% Setting up the matrix describing the active Piezoelectric layer

% The material parameters

rho_pz = 7800;  
c_pz = 4613;    
l_pz = 2*10^-3;
r_pz = 2.5*10^-3;
A_pz = pi*r_pz^2;
Z_0_pz = rho_pz * c_pz * A_pz;
f_0_pz = c_pz/(2*l_pz);
omg_0_pz = pi*c_pz/l_pz;
gamma_pz = pi*f/f_0_pz;
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
C_0_pz = A_pz*eps_33/l_pz;
k_pz = 0.46;                        % coupling factor in thickness direction
phi = k_pz *sqrt(Z_0_pz*C_0_pz*omg_0_pz/pi);


z_b = Z_b/Z_0_pz;
Q = cos(gamma_pz) - 1 + j*z_b*sin(gamma_pz);
pz_frac = 1/(phi*Q);
pz_mat_1 = [1 (j*phi^2)/(omg*C_0_pz); ...           % Revisit as the fraction of entry 2,1 is not clear!!
            j*omg*C_0_pz 0];
pz_mat_2 = [cos(gamma_pz)+j*z_b*sin(gamma_pz) Z_0_pz*(z_b*cos(gamma_pz)+j*sin(gamma_pz)); ...
            (j*sin(gamma_pz))/Z_0_pz 2*(cos(gamma_pz)-1)+j*z_b*sin(gamma_pz)];     
pz_mat = pz_frac*pz_mat_1*pz_mat_2;
%% Setting up the final transducer representing matrix

% Here the piezoelectric matric is multiplied by the matrix obtained from
% the transmission intermediate layers:

Final_trans_mat = pz_mat * Trans_mat;

%% Defining the sensitivity functions
A = Final_trans_mat(1,1);
B = Final_trans_mat(1,2);
C = Final_trans_mat(2,1);
D = Final_trans_mat(2,2);

Z_FP = Z_0_FeC;
Z_E = (A*Z_FP + B)/(C*Z_FP + D);            % electrical input impedance
V_IL = Z_FP/(A*Z_FP + B);                   % voltage transfer ratio between input voltage and output force


%% Plotting the resulting pressure field

V_M = 150;                                                                          % Input voltage
c_fluid = 1300;
lambda =c_fluid/f;                                                                  % wavelength of the pressure field
F_surface = V_M * V_IL;                                                       % Total force on trasnducer surface
P_surface = F_surface/A_FP;                                                  % pressure on the transducer surface

Kin_Vis=32*10^-6;
damp_coeff =(2*Kin_Vis*omg^2)/(3*c_fluid^3);                                          % Damping coefficient for attenuation within the liquid
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

t = linspace(0,10*pi/omg,1000);                                               % time duration during which the pressure field
                                                                            % is monitored.
x = linspace(0,5*0.0013,1000);                                              % Distance over which the pressure field is monitored.

% Here the P_field matrix is filled with the pressure amplitudes and phase
% shifts at a certain location in space and moment in time.
 for n = 1:1000
     for m = 1:1000
     P_field(n,m) = P_surface*(R1+1) * exp(1i*omg*t(m)-1i*x(n)*2*pi/lambda) ...
         - R2*P_surface * exp(1i*omg*t(m)+1i*(x(n))*2*pi/lambda);
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


