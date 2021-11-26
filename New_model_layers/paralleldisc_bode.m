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

f = linspace(1,2*10^6,100000);      
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
Back_mat_Cu = zeros(2*length(f),2);

for i = 1:length(f)
row1 = 2*i-1;
row2 = 2*i;
Back_mat_Cu(row1:row2,1:2) = [cos(gamma_Cu(i)) j*Z_0_Cu*sin(gamma_Cu(i)); ...
               (j*sin(gamma_Cu(i)))/Z_0_Cu cos(gamma_Cu(i))];
end
          
           
           
% Similarily the matrix is set up for the thick, cone shaped backing
% material

% the material parameters

rho_FeC = 8000;  
c_FeC = 5800;    
r_FeC = 2.5*10^-3;
A_FeC = pi*r_FeC^2;
Z_0_FeC = rho_FeC * c_FeC * A_FeC;
f_0_FeC = c_FeC/(2);                     % As the thickness of the backing is assumed to approach infinity, it is omitted. Not sure about this! Check!!
gamma_FeC = pi*f/f_0_FeC;

% the matrix
Back_mat_FeC = zeros(2*length(f),2);

for j = 1:length(f)
row1 = 2*j-1;
row2 = 2*j;
  
Back_mat_FeC(row1:row2,1:2) = [cos(gamma_FeC(j)) j*Z_0_FeC*sin(gamma_FeC(j)); ...
                (j*sin(gamma_FeC(j)))/Z_0_FeC cos(gamma_FeC(j))];
end

rho_glue = 960;  
c_glue = 1220;    % unsure about this!!
l_glue = 0.06*10^-6;    % to be varied.
r_glue = 2.5*10^-3;
A_glue = pi*r_Cu^2;
Z_0_glue = rho_Cu * c_Cu * A_Cu;
f_0_glue = c_Cu/(2*l_Cu);
gamma_glue = pi*f/f_0_Cu;

glue_mat = zeros(2*length(f),2);

for han = 1:length(f)
row1 = 2*han-1;
row2 = 2*han;

glue_mat(row1:row2,1:2) = [cos(gamma_glue(han)) j*Z_0_glue*sin(gamma_glue(han)); ...
               (j*sin(gamma_glue(han)))/Z_0_glue cos(gamma_glue(han))];
end
          
% The matrices are multiplied and the total impedance Z_b of the backing
% layers is obtained.
Back_mat = zeros(2*length(f),2);
Z_b = zeros(length(f),1);

for ii = 1:length(f)
row1 = 2*ii-1;
row2 = 2*ii;
    
Back_mat(row1:row2,1:2) = glue_mat(row1:row2,1:2)*Back_mat_Cu(row1:row2,1:2) ...
    *glue_mat(row1:row2,1:2)*Back_mat_FeC(row1:row2,1:2);
Back_mat_current = Back_mat(row1:row2,1:2);
A_b = Back_mat_current(1,1);
B_b = Back_mat_current(1,2);
C_b = Back_mat_current(2,1);
D_b = Back_mat_current(2,2);

Z_0b = Z_0_FeC;
Z_b(ii) = (A_b*Z_0b + B_b)/(C_b*Z_0b + D_b);    % Here Z_0b is the impedance of the final backing layer only!
end

%% Setting up the transmission layer matrix

% The intermediate layers in the transmission direction are obtained in the
% same way. Again starting with the intermediate layer closest to the
% peizoelectric active layer.

% The material parameters were already given before!
Trans_mat_Cu = zeros(length(f),2);

for jj = 1:length(f)
row1 = 2*jj-1;
row2 = 2*jj;
   
Trans_mat_Cu(row1:row2,1:2) = [cos(gamma_Cu(jj)) j*Z_0_Cu*sin(gamma_Cu(jj)); ...
               (j*sin(gamma_Cu(jj)))/Z_0_Cu cos(gamma_Cu(jj))];
end

% For now the glue layer is omitted and the matrix is setup for the front
% plate only. Which is made out of Steel, just as the backing material.

% material properties of steel front plate.
r_FP = 4.5*10^-3;
A_FP = pi*r_FP^2;
Z_0_FP = rho_FeC * c_FeC * A_FP;
l_FeC = 0.5*10^-3;
f_0_FP = c_FeC/(2*l_FeC);                    
gamma_FeC = pi*f/f_0_FP;

Trans_mat_FP = zeros(2*length(f),2);

for x = 1:length(f)
row1 = 2*x-1;
row2 = 2*x;

Trans_mat_FP(row1:row2,1:2) = [cos(gamma_FeC(x)) j*Z_0_FP*sin(gamma_FeC(x)); ...
                (j*sin(gamma_FeC(x)))/Z_0_FP cos(gamma_FeC(x))];
end

%The housing, simply modelled as a another transmission plate

rho_house = rho_FeC;  
c_house = c_FeC;    
l_house = 3.25*10^-3;    
r_house = 3.5*10^-3;
A_house = pi*r_Cu^2;
Z_0_house = rho_Cu * c_Cu * A_Cu;
f_0_house = c_Cu/(2*l_Cu);
gamma_house = pi*f/f_0_Cu;

house_mat = zeros(2*length(f),2);

for kla = 1:length(f)
    row1 = 2*kla-1;
    row2 = 2*kla;
    
house_mat(row1:row2,1:2) = [cos(gamma_house(kla)) j*Z_0_house*sin(gamma_house(kla)); ...
               (j*sin(gamma_house(kla)))/Z_0_house cos(gamma_house(kla))];
end
% Multiplication of the transmission piezoelectric inactive layers gives:
Trans_mat = zeros(2*length(f),2);

for xx = 1:length(f)
row1 = 2*xx-1;
row2 = 2*xx;

Trans_mat(row1:row2,1:2) = glue_mat(row1:row2,1:2)*Trans_mat_Cu(row1:row2,1:2) ...
    * glue_mat(row1:row2,1:2)*Trans_mat_FP(row1:row2,1:2)*house_mat(row1:row2,1:2);
end
                   
%% Setting up the matrix describing the active Piezoelectric layer

% The material parameters

rho_pz = 7800;  
c_pz = 4613;    
l_pz = 2*10^-3;
r_pz = 2.5*10^-3;
A_pz = pi*r_pz^2;
Z_0_pz = rho_pz * c_pz * A_pz;
N = 2;
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
C_s = A_pz*eps_33/l_pz;
s_33 = 14.2*10^-12;
d_33 = 265*10^-12;
h_33 = d_33/(s_33*eps_33);
h = h_33;
Z = A_pz*rho_pz*c_pz;
Z_inv = 1/Z;
theta = omg*l_pz/c_pz;
k = 0.46;
sigma = k^2./theta;
phi = acos((cos(theta)-sigma.*sin(theta))./(1-sigma.*sin(theta)));
R = (sqrt(sin(theta)-2*sigma.*(1-cos(theta))))./sin(theta);
R_inv = 1./R;

T_11 = cos(N*phi);
T_12 = -j*Z*R.*sin(N*phi);
T_13 = -h*C_s*tan(0.5*phi).*sin(N*phi);
T_14 = 0;
T_21 = -j*Z_inv.*R_inv.*sin(N*phi);
T_22 = cos(N*phi);
T_23 = -j*h*C_s.*Z_inv.*R_inv.*tan(0.5*phi).*(cos(N*phi)-(-1)^N);
T_24 = 0;
T_31 = 0;
T_32 = 0;
T_33 = (-1)^N;
T_34 = 0;
T_41 = -j*h*C_s.*Z_inv.*R_inv.*tan(0.5*phi).*(cos(N*phi)-(-1)^N);
T_42 = -h*C_s*tan(0.5*phi).*sin(N*phi);
T_43 = j*((N*(-1)^N)*(1+2*sigma.*R_inv.*tan(0.5*phi))+sigma.*R_inv.*((tan(0.5*phi)).^2).*sin(N*phi)).*omg*C_s;
T_44 = (-1)^N;

z_b = Z_b/Z_0_pz;

A_c = zeros(length(f),1);
B_c = zeros(length(f),1);
C_c = zeros(length(f),1);
D_c = zeros(length(f),1);
pz_mat = zeros(2*length(f),2);

for y = 1:length(f)
row1 = 2*y-1;
row2 = 2*y;

A_c(y) = T_31 - T_33*(T_21(y)*z_b(y) + T_11(y))/(T_23(y)*z_b(y) + T_13(y));
B_c(y) = T_32 - T_33*(T_22(y)*z_b(y) + T_12(y))/(T_23(y)*z_b(y) + T_13(y));
C_c(y) = T_41(y) - T_43(y)*(T_21(y)*z_b(y) + T_11(y))/(T_23(y)*z_b(y) + T_13(y));
D_c(y) = T_42(y) - T_43(y)*(T_22(y)*z_b(y) + T_12(y))/(T_23(y)*z_b(y) + T_13(y));
pz_mat(row1:row2,1:2) = [A_c(y) B_c(y); C_c(y) D_c(y)];
end
%% Setting up the final transducer representing matrix

% Here the piezoelectric matric is multiplied by the matrix obtained from
% the transmission intermediate layers:

Final_trans_mat = zeros(2*length(f),2);

for yy = 1:length(f)
row1 = 2*yy-1;
row2 = 2*yy;

Final_trans_mat(row1:row2,1:2) = pz_mat(row1:row2,1:2) * Trans_mat(row1:row2,1:2);
end

%% Defining the sensitivity functions
Z_E = zeros(length(f),1);
V_IL = zeros(length(f),1);

for g = 1:length(f)
row1 = 2*g-1;
row2 = 2*g;
Final_trans_mat1 = Final_trans_mat(row1:row2,1:2);
A = Final_trans_mat1(1,1);
B = Final_trans_mat1(1,2);
C = Final_trans_mat1(2,1);
D = Final_trans_mat1(2,2);

Z_FP = Z_0_FeC;
Z_E(g) = (A*Z_FP + B)/(C*Z_FP + D);            % electrical input impedance
V_IL(g) = Z_FP/(A*Z_FP + B);                   % voltage transfer ratio between input voltage and output force
end

admittance = 1./Z_E;
susceptance = imag(admittance);
phase_admittance = rad2deg(angle(admittance));
mag_admittance = 20*log10(abs(admittance));

phase_impedance = rad2deg(angle(Z_E));
phase_tf = rad2deg(angle(V_IL));
mag_impedance = 20*log10(abs(Z_E));
mag_tf = 20*log10(abs(V_IL));

%[pks1 ind1] = findpeaks(mag_impedance,'MinPeakDistance',5000,'MinPeakProminence',0.01);
%plot(f(ind1),pks1,'or'); hold on; plot(f(ind),pks,'or'); hold on;
figure
subplot(2,1,1); plot(f,mag_impedance); title('Bode magnitude plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]'); 
subplot(2,1,2); plot(f,phase_impedance); title('Phase plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');


% [pks ind] = findpeaks(mag_tf,'MinPeakDistance',5000,'MinPeakProminence',0.01);

figure
subplot(2,1,1);  plot(f,mag_tf); title('Bode magnitude plot of the transfer function V_{in} --> F_{out}'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');
subplot(2,1,2); plot(f,phase_tf); title('Phase plot of the transfer function V_{in} --> F_{out}'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');

figure
subplot(2,1,1); plot(f,mag_admittance); title('Bode magnitude plot of the admittance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');
subplot(2,1,2); plot(f,phase_admittance); title('Phase plot of the admittance'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');

figure
plot(f,susceptance); title('Bode magnitude plot of the susceptance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');

