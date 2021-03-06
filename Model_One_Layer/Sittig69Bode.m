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

f = linspace(1,2*10^6,1000); 

omg = 2*pi*f;
j = sqrt(-1);
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

for ja = 1:length(f)
row1 = 2*ja-1;
row2 = 2*ja;
Back_mat_Cu(row1:row2,1:2) = [cos(gamma_Cu(ja)) j*Z_0_Cu*sin(gamma_Cu(ja)); ...
               (j*sin(gamma_Cu(ja)))/Z_0_Cu cos(gamma_Cu(ja))];
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

for jj = 1:length(f)
row1 = 2*jj-1;
row2 = 2*jj;
  
Back_mat_FeC(row1:row2,1:2) = [cos(gamma_FeC(jj)) j*Z_0_FeC*sin(gamma_FeC(jj)); ...
                (j*sin(gamma_FeC(jj)))/Z_0_FeC cos(gamma_FeC(jj))];
end

rho_glue = 960;  
c_glue = 1220;    % unsure about this!!
l_glue = 0.06*10^-6;    % to be varied.
r_glue = 2.5*10^-3;
A_glue = pi*r_glue^2;
Z_0_glue = rho_glue * c_glue * A_glue;
f_0_glue = c_glue/(2*l_glue);
gamma_glue = pi*f/f_0_glue;

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
A_house = pi*r_house^2;
Z_0_house = rho_house * c_house * A_house;
f_0_house = c_house/(2*l_house);
gamma_house = pi*f/f_0_house;

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
f_0_pz = c_pz/(2*l_pz);
omg_0_pz = pi*c_pz/l_pz;
gamma_pz = pi*f/f_0_pz;
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
C_0_pz = A_pz*eps_33/l_pz;
k_pz = 0.46;                        % coupling factor in thickness direction
phi = k_pz *sqrt(Z_0_pz*C_0_pz*omg_0_pz/pi);


z_b = Z_b/Z_0_pz;
Q1 = (cos(gamma_pz) - 1);
Q2 =  j*z_b(:);
Q3 = Q2.*sin(gamma_pz(:));
Q = Q1(:)+Q3;
pz_frac = 1./(phi*Q);

pz_mat_1_12 = (j*phi^2)./(omg*C_0_pz);
pz_mat_1_21 = j*omg*C_0_pz;




pz_mat_1 = zeros(2*length(f),2);
pz_mat_2 = zeros(2*length(f),2);
pz_mat = zeros(2*length(f),2);

for y = 1:length(f)
row1 = 2*y-1;
row2 = 2*y;
pz_mat_1(row1:row2,1:2) = [1 pz_mat_1_12(y); ...           % Revisit as the fraction of entry 2,1 is not clear!!
            pz_mat_1_21(y) 0];
pz_mat_1_current = pz_mat_1(row1:row2,1:2);

pz_mat_2_11 = cos(gamma_pz(y))+j*z_b(y).*sin(gamma_pz(y));
pz_mat_2_12 = Z_0_pz*(z_b(y).*cos(gamma_pz(y))+j*sin(gamma_pz(y)));
pz_mat_2_21 = (j*sin(gamma_pz(y)))/Z_0_pz;
pz_mat_2_22 = 2*(cos(gamma_pz(y))-1)+j*z_b(y).*sin(gamma_pz(y));
pz_mat_2(row1:row2,1:2) = [pz_mat_2_11 pz_mat_2_12; ...
            pz_mat_2_21 pz_mat_2_22];
pz_mat_2_current = pz_mat_2(row1:row2,1:2);
        
pz_mat(row1:row2,1:2) = pz_frac(y)*pz_mat_1_current*pz_mat_2_current;
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
k_c = omg/(2/3*3*10^8);                                %Wave number in cable
l_c = 1;                                              %Length of cable
Z_c = j*4.3*10^-6*omg;                                      %Impedance of cable
T_c_11 = cos(k_c*l_c);                                   %Transfer matrix of cable
T_c_12 = -j*Z_c.*sin(k_c*l_c);                         %Transfer matrix of cable

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
%Z_E(g) = (A*Z_FP + B)/(C*Z_FP + D);            % electrical input impedance
Z_E(g) = ((A*Z_FP + B)*T_c_11(g)+T_c_12(g))/(C*Z_FP + D);
%V_IL(g) = Z_FP/(A*Z_FP + B);                   % voltage transfer ratio between input voltage and output force
V_IL(g) = Z_FP/((A*Z_FP + B)*T_c_11(g)+T_c_12(g));
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
plot(f,imag(Z_E)); title('Bode magnitude plot of the susceptance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]'); ylim([-0.2 0.2]*10^6)

