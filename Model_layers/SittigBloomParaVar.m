% This script is an intent to accurately model a piezoelectric transducer,
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

% Be careful when changing the material parameters!! For now some
% parameters can be changed, however, the order is very important! Per
% changeable parameter, there are 10 different possibilities, resulting in
% 10*1000 = 10000 loops. In the first loop the frequency response for ONE
% value of the parameter is used. This counts for every variable material
% parameter! Therefore, if e.g. both l_glue and l_Cu are changed to be of
% the form linspace(x,y,10), this would mean that for 1000 loops the first
% entries of both vectors are utilized and then consistantly the next
% entries are used for the coming 1000 loops. Hence, if one would like to
% see the effect of a change in one parameter only, the values of all the
% other parameters have to be set to one value, in the form:
% ones(10,1)*specific_value.

%% close all former work
clear all; close all; clc;

%% Parameter initialization

f = linspace(1,2*10^6,1000);      
omg = 2*pi*f;
j = 1i;
%% Setting up the backing material matrix

% First the layer closest to the piezoelectric layer, hence the copper
% electrode.

% the material parameters

rho_Cu = 8960;  % Pure! copper density.
c_Cu = 3570;    % pure copper speed of sound.

% depending on the type of parameter variation, one of the following
% options for l_Cu should be chosen.
l_Cu = ones(10,1)*0.5e-3;   % choose this if the effect of one other parameter variation is to be obtained.                       
% l_Cu = linspace(0.4,0.6,10)*1e-3;   % Choose this to either see the effect of variations in copper thickness, or to combine parameter variation.
% l_Cu = l_Cu(:);
r_Cu = 2.5*10^-3;
A_Cu = pi*r_Cu^2;
Z_0_Cu = rho_Cu * c_Cu * A_Cu;
f_0_Cu = c_Cu./(2*l_Cu);
gamma_Cu = pi*f./f_0_Cu;

% The matrix
Back_mat_Cu = zeros(2*length(f),2*length(l_Cu));

for ta = 1:length(l_Cu)
    col1_Cu = 2*ta-1;
    col2_Cu = 2*ta;
for i = 1:length(f)
row1 = 2*i-1;
row2 = 2*i;
Back_mat_Cu(row1:row2,col1_Cu:col2_Cu) = [cos(gamma_Cu(ta,i)) j*Z_0_Cu*sin(gamma_Cu(ta,i)); ...
               (j*sin(gamma_Cu(ta,i)))/Z_0_Cu cos(gamma_Cu(ta,i))];
end
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

for is = 1:length(f)
row1 = 2*is-1;
row2 = 2*is;
  
Back_mat_FeC(row1:row2,1:2) = [cos(gamma_FeC(is)) j*Z_0_FeC*sin(gamma_FeC(is)); ...
                (j*sin(gamma_FeC(is)))/Z_0_FeC cos(gamma_FeC(is))];
end

rho_glue = 960;  
c_glue = 1220;    % unsure about this!!
l_glue = ones(10,1)*100e-6;    % choose this if the effect of one other parameter variation is to be obtained.
% l_glue = linspace(10,100,10)*10^-6; % Choose this to either see the effect of variations in glue thickness, or to combine parameter variation.
% l_glue = l_glue(:);
r_glue = 2.5*10^-3;
A_glue = pi*r_glue^2;
Z_0_glue = rho_glue * c_glue * A_glue;
f_0_glue = c_glue./(2*l_glue);

% material properties of steel front plate.
r_FP = 4.5*10^-3;
A_FP = pi*r_FP^2;
Z_0_FP = rho_FeC * c_FeC * A_FP;
l_FeC = 0.5*10^-3;
f_0_FP = c_FeC/(2*l_FeC);                    
gamma_FeC = pi*f/f_0_FP;

rho_house = rho_FeC;  
c_house = c_FeC;    
l_house = ones(10,1)*3.25*10^-3;   % choose this if the effect of one other parameter variation is to be obtained. 
% l_house = linspace(0.5,1.5,10)*3.25*10^-3;  % Choose this to either see the effect of variations in housing thickness, or to combine parameter variation.
% l_house = l_house(:);
r_house = 3.5*10^-3;
A_house = pi*r_house^2;
Z_0_house = rho_house * c_house * A_house;
f_0_house = c_house./(2*l_house);
gamma_house = pi*f./f_0_house;

% pz para
rho_pz = 7800;  
c_pz = 4613;    
l_pz = ones(10,1)*2e-3;  % choose this if the effect of one other parameter variation is to be obtained.
% l_pz = linspace(0.5,1.5,10)*2e-3; % Choose this to either see the effect of variations in piezo thickness, or to combine parameter variation.
% l_pz = l_pz(:);
r_pz = 2.5*10^-3;
A_pz = pi*r_pz^2;
Z_0_pz = rho_pz * c_pz * A_pz;
N = 2;
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
C_s = A_pz*eps_33./l_pz;
s_33 = 14.2*10^-12;
d_33 = 265*10^-12;
h_33 = d_33/(s_33*eps_33);
h = h_33;
Z = A_pz*rho_pz*c_pz;
Z_inv = 1/Z;
theta = omg.*l_pz./c_pz;
k = 0.46;
sigma = k^2./theta;
phi = acos((cos(theta)-sigma.*sin(theta))./(1-sigma.*sin(theta)));
R = (sqrt(sin(theta)-2*sigma.*(1-cos(theta))))./sin(theta);
R_inv = 1./R;
T_11 = cos(N*phi);
T_12 = -j*Z*R.*sin(N*phi);
T_13 = -h*C_s.*tan(0.5*phi).*sin(N*phi);
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
T_42 = -h*C_s.*tan(0.5*phi).*sin(N*phi);
T_43 = j*((N*(-1)^N)*(1+2*sigma.*R_inv.*tan(0.5*phi))+sigma.*R_inv.*((tan(0.5*phi)).^2).*sin(N*phi)).*omg.*C_s;
T_44 = (-1)^N;
inductance_cable = 4.3*10^-6;

glue_mat = zeros(2*length(f),2*length(l_glue));
Back_mat = zeros(2*length(f),2*length(l_glue));
Z_b = zeros(length(f),length(l_glue));
Trans_mat = zeros(2*length(f),2*length(l_glue));
Final_trans_mat = zeros(2*length(f),2*length(l_glue));
Z_E = zeros(length(f),length(l_glue));
V_IL = zeros(length(f),length(l_glue));
phase_tf = zeros(length(f),length(l_glue));

for er = 1:length(l_glue)
col1 = 2*er-1;
col2 = 2*er;
gamma_glue = pi*f/f_0_glue(er);
for han = 1:length(f)
row1 = 2*han-1;
row2 = 2*han;
    glue_mat(row1:row2,col1:col2) = [cos(gamma_glue(han)) j*Z_0_glue*sin(gamma_glue(han)); ...
               (j*sin(gamma_glue(han)))/Z_0_glue cos(gamma_glue(han))];
           
% The matrices are multiplied and the total impedance Z_b of the backing
% layers is obtained.
Back_mat(row1:row2,col1:col2) = glue_mat(row1:row2,col1:col2)*Back_mat_Cu(row1:row2,col1:col2) ...
    *glue_mat(row1:row2,col1:col2)*Back_mat_FeC(row1:row2,1:2);
Back_mat_current = Back_mat(row1:row2,col1:col2);
A_b = Back_mat_current(1,1);
B_b = Back_mat_current(1,2);
C_b = Back_mat_current(2,1);
D_b = Back_mat_current(2,2);

Z_0b = Z_0_FeC;
Z_b(han,er) = (A_b*Z_0b + B_b)/(C_b*Z_0b + D_b);    % Here Z_0b is the impedance of the final backing layer only!

Trans_mat_Cu = Back_mat_Cu(row1:row2,col1:col2);
Trans_mat_FP = [cos(gamma_FeC(han)) j*Z_0_FP*sin(gamma_FeC(han)); ...
                (j*sin(gamma_FeC(han)))/Z_0_FP cos(gamma_FeC(han))];
house_mat = [cos(gamma_house(er,han)) j*Z_0_house*sin(gamma_house(er,han)); ...
               (j*sin(gamma_house(er,han)))/Z_0_house cos(gamma_house(er,han))];
Trans_mat(row1:row2,col1:col2) = glue_mat(row1:row2,col1:col2)*Trans_mat_Cu ...
    * glue_mat(row1:row2,col1:col2)*Trans_mat_FP*house_mat;
z_b = Z_b(han,er)/Z_0_pz;

A_c = T_31 - T_33*(T_21(er,han)*z_b + T_11(er,han))/(T_23(er,han)*z_b + T_13(er,han));
B_c = T_32 - T_33*(T_22(er,han)*z_b + T_12(er,han))/(T_23(er,han)*z_b + T_13(er,han));
C_c = T_41(er,han) - T_43(er,han)*(T_21(er,han)*z_b + T_11(er,han))/(T_23(er,han)*z_b + T_13(er,han));
D_c = T_42(er,han) - T_43(er,han)*(T_22(er,han)*z_b + T_12(er,han))/(T_23(er,han)*z_b + T_13(er,han));
pz_mat = [A_c B_c; C_c D_c];

Final_trans_mat(row1:row2,col1:col2) = pz_mat*Trans_mat(row1:row2,col1:col2);
Final_trans_mat1 = Final_trans_mat(row1:row2,col1:col2);
A = Final_trans_mat1(1,1);
B = Final_trans_mat1(1,2);
C = Final_trans_mat1(2,1);
D = Final_trans_mat1(2,2);
Z_FP = Z_0_FeC;
Z_E(han,er) = -(A*Z_FP + B)/(C*Z_FP + D);            % electrical input impedance
V_IL(han,er) = Z_FP/(A*Z_FP + B);   
phase_tf(han,er) = rad2deg(phase(V_IL(han,er)));
end
end         

admittance = 1./Z_E;
susceptance = imag(admittance);
conductance = real(admittance);
phase_admittance = -rad2deg(angle(admittance));
mag_admittance = 20*log10(abs(admittance));
mag_conductance = 20*log10(abs(conductance));
phase_impedance = -rad2deg(angle(Z_E));



% first_quadrant = find(real(V_IL)>0 & imag(V_IL)>0);
% second_quadrant = find(real(V_IL)<0 & imag(V_IL)>0);
% third_quadrant = find(real(V_IL)<0 & imag(V_IL)<0);
% fourth_quadrant = find(real(V_IL)>0 & imag(V_IL)<0);
% quadrant_order = [first_quadrant; second_quadrant; third_quadrant; fourth_quadrant];
% phase1 = rad2deg(angle(V_IL(first_quadrant)));
% phase2 = rad2deg(angle(V_IL(second_quadrant)));
% phase3 = 360 + rad2deg(angle(V_IL(third_quadrant)));
% phase4 = 360 + rad2deg(angle(V_IL(fourth_quadrant)));
% phasers = [phase1; phase2; phase3; phase4];
% sync_vec = [quadrant_order phasers];
% [a_sorted, a_order] = sort(sync_vec(:,1));
% B = sync_vec(:,2);
% newB = B(a_order,:);
% sortedvector = [a_sorted newB];


phase_resistance = rad2deg(angle(real(Z_E)));
mag_impedance = 20*log10(abs(Z_E));
mag_tf = 20*log10(abs(V_IL));
mag_resistance = 20*log10(abs(real(Z_E)));


ds = tabularTextDatastore('C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\newData\project','FileExtensions','.csv');
T = readall(ds);
ind_freq = 13:94:93919;
ind_phase = 27:94:93933;
ind_ptp1 = 55:94:93961;
ind_ptp2 = 69:94:93975;
ind_ptp4 = 83:94:93989;
frequencies = T(ind_freq,1);
phase_shift = T(ind_phase,1);
ptp_voltage = T(ind_ptp1,1);
ptp_current = T(ind_ptp2,1);
ptp_trans2 = T(ind_ptp4,1);

freq = table2array(frequencies);
freq = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),freq);
phase_shift = table2array(phase_shift);
phase_shift = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),phase_shift);
ptp_voltage = table2array(ptp_voltage);
ptp_voltage = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_voltage);
ptp_current = table2array(ptp_current);
ptp_current = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_current);
ptp_trans2 = table2array(ptp_trans2);
ptp_trans2 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_trans2);

es = tabularTextDatastore('C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\newData\projectcontinued');
S = readall(es);

ind_freq1 = 13:94:37613;
ind_phase1 = 27:94:37627;
ind_ptp1 = 55:94:37655;
ind_ptp2 = 69:94:37669;
ind_ptp4 = 83:94:37683;
frequencies1 = S(ind_freq1,1);
phase_shift1 = S(ind_phase1,1);
ptp_voltage1 = S(ind_ptp1,1);
ptp_current1 = S(ind_ptp2,1);
ptp_trans21 =  S(ind_ptp4,1);

freq1 = table2array(frequencies1);
freq1 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),freq1);
phase_shift1 = table2array(phase_shift1);
phase_shift1 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),phase_shift1);
ptp_voltage1 = table2array(ptp_voltage1);
ptp_voltage1 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_voltage1);
ptp_current1 = table2array(ptp_current1);
ptp_current1 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_current1);
ptp_trans21 = table2array(ptp_trans21);
ptp_trans21 = cellfun(@(x)str2double(regexp(x,'\d*\.\d*','match')),ptp_trans21);

freq_vec = [freq; freq1];
multiplyer = freq_vec(1:901)*1000;
multiplyer2 = freq_vec(902:end)*1000000;
freq_vec = [multiplyer; multiplyer2];
phase_shift_vec = [phase_shift; phase_shift1];
ptp_voltage_vec = [ptp_voltage; ptp_voltage1];
ptp_current_vec = [ptp_current; ptp_current1];
ptp_trans2_vec = [ptp_trans2; ptp_trans21];

j = 1i;
current_phasor_amp = ptp_current_vec./2000;             % amplitude in mA.
current_angle = deg2rad(phase_shift_vec);
current_phasor = current_phasor_amp.*exp(j*current_angle);

impedance_phasor = (ptp_voltage_vec./2)./current_phasor;

phase_impedance2 = rad2deg(angle(impedance_phasor));
% high_phase = find(phase_impedance2(phase_impedance2>=90));
% phase_impedance2(high_phase) = 180-phase_impedance2(high_phase);
phase_resistance2 = rad2deg(angle(real(impedance_phasor)));
mag_impedance2 = 20*log10(abs(impedance_phasor));
mag_resistance2 = 20*log10(abs(real(impedance_phasor)));


admittance2 = 1./impedance_phasor;
phase_admittance2 = rad2deg(angle(admittance2));
mag_admittance2 = 20*log10(abs(admittance2));
susceptance2 = imag(admittance2);
conductance2 = real(admittance2);
mag_conductance2 = 20*log10(abs(conductance2));

% [pks1 ind1] = findpeaks(mag_impedance,'MinPeakDistance',5000,'MinPeakProminence',0.01);
% plot(f(ind1),pks1,'or'); hold on; plot(f(ind),pks,'or'); hold on;
figure(1)
hold on
figure(2)
hold on
figure(3)
hold on
figure(4)
hold on
figure(5)
hold on
figure(6)
hold on

for uu = 1:length(l_glue)
% row_begin = length(f)*uu-(length(f)-1);
% row_end = length(f)*uu;

figure(1);
subplot(2,1,1); plot(f,mag_impedance(:,uu)); hold on; plot(freq_vec,mag_impedance2,'r'); title('Frequency response of the impedance'); ...
xlabel('frequency [Hz]'); ylabel('magnitude [DB]'); legend('Theoretical system impedance'); 
conversion_vec = 180*ones(length(freq_vec),1);
subplot(2,1,2); plot(f,phase_impedance(:,uu)); hold on; plot(freq_vec,phase_impedance2-conversion_vec,'r');  ...
xlabel('frequency [Hz]'); ylabel('phase [degree]'); legend('Theoretical system impedance phase'); ylim([-180 180]);

figure(2);
subplot(2,1,1);  plot(f./1e6,mag_tf(:,uu)); hold on; title('Frequency response of S_{FV}'); ...
    xlabel('frequency [MHz]'); ylabel('magnitude [DB]'); legend('Theoretical system frequency response'); ylim([-250 50])
subplot(2,1,2); plot(f./1e6,phase_tf(:,uu)); hold on; ...
    xlabel('frequency [MHz]'); ylabel('phase [degree]'); legend('Theoretical system frequency response'); 

figure(3);
subplot(2,1,1); plot(f./1e6,mag_admittance(:,uu)); hold on; plot(freq_vec./1e6,mag_admittance2,'r');title('Bode magnitude plot of the admittance'); ...
    xlabel('frequency [MHz]'); ylabel('magnitude [DB]'); legend('Theoretical system admittance');
subplot(2,1,2); plot(f./1e6,phase_admittance(:,uu));  hold on; plot(freq_vec./1e6,-phase_admittance2,'r');
 title('Theoretical system admittance phase'); ...
    xlabel('frequency [MHz]'); ylabel('phase [degree]'); ylim([-150 150])

figure(4)
plot(f,susceptance(:,uu)); hold on; plot(freq_vec,susceptance2,'r'); title('Bode magnitude plot of the susceptance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');

figure(5)
semilogy(f,mag_conductance(:,uu)); hold on; semilogy(freq_vec,mag_conductance2); title('magnitude plot of the conductance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');

figure(6)
semilogy(f,mag_resistance(:,uu)); hold on; semilogy(freq_vec,mag_resistance2); title('magnitude plot of the resistance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');

 
end







