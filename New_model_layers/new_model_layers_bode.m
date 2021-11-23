close all;
clear all;
clc;
%% Parameter initialization

j = 1i;                                                                     % imaginary value
f = linspace(1,2*10^6,10000);                                                               % operating frequency [Hz]
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
k_n_squared = h_33^2 *eps_33/c_33 * ones(length(omega),1);
k_n_squared = k_n_squared';
k_squared = k_n_squared;
gamma_sinus = sin(gamma);
s_intermediate = k_squared.*gamma_sinus;
s = s_intermediate./gamma;
%gamma_mn = [pi*f/f_mn pi*f/f_backing];
r = 2.5*10^-3; 
S = pi*r^2;
C_0= S*eps_33/(d_n);
Z_0 = rho_n*v_n*S;
omega_mn = 2*pi*f_mn;
Z_mn = Z_0;
phi = sqrt(k_squared)*sqrt(omega_mn*c_mn*Z_mn/pi);
gamma_cosinus = ones(1,length(gamma)) - cos(gamma);
c_intermediate = k_squared .* gamma_cosinus;
c = c_intermediate./gamma;
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




%% 4 port model matrix


T_11 = (cos(gamma) - s)./(1-s);
T_12 = (j*Z_0*(sin(gamma) - 2*c))./(1-s);
cosinus_gamma = cos(gamma)-1;
T_13 = -cosinus_gamma.*phi./(1-s);
T_14 = 0;
T_21 = (j*sin(gamma))./(Z_0*(1-s));
T_22 = (cos(gamma) - s)./(1-s);
T_23 = -(j*phi.*sin(gamma))./(1-s);
T_24 = 0;
T_31 = 0;
T_32 = 0;
T_33 = 1;
T_34 = 0;
T_41_intermediate = -(j*sin(gamma))./(Z_0*(1-s));
T_41 = T_41_intermediate .* phi;
T_42 = -((cos(gamma) - 1).*phi)./(1-s);
T_43 = (j*omega*C_0)./(1-s);
T_44 = 1;

%% Intermediate layers (ex. glue layers, copper electrode, etceteraness)

r_copper = 2.5*10^-3;
A_copper = pi*r_copper^2;
rho_copper = 8920;
c_copper = 3570;
d_copper = 0.5*10^-3;
v_copper = sqrt(c_copper/rho_copper);

Z_copper = A_copper*rho_copper*c_copper;
gamma_copper = omega * d_copper/v_copper;

A_copper = zeros(2*length(omega),2);

for i = 1:length(omega)
    entry1 = 2*i - 1;
    entry2 = 2*i;
A_copper(entry1:entry2,1:2) = [cos(gamma_copper(i)) j*Z_copper*sin(gamma_copper(i)); ...
    (j*sin(gamma_copper(i)))/Z_copper cos(gamma_copper(i))];
end

rho_FP = 7850;
r_FP = 4.5*10^-3;
c_FP = 3230;
S_FP = pi*r_FP^2;
d_FP = 0.5*10^-3;
v_FP = sqrt(c_FP/rho_FP);

Z_FP = rho_FP*c_FP*S_FP;
gamma_FP = omega*d_FP/v_FP;

A_FP = zeros(2*length(omega),2);

for ii = 1:length(omega)
    entry1 = 2*ii - 1;
    entry2 = 2*ii;
A_FP(entry1:entry2,1:2) = [cos(gamma_FP(ii)) j*Z_FP*sin(gamma_FP(ii)); ...
    (j*sin(gamma_FP(ii)))/Z_FP cos(gamma_FP(ii))];
end

Z_b = zeros(length(omega),1);

for j = 1:length(omega)
    entry1 = 2*j - 1;
    entry2 = 2*j;
A_m = A_copper;
A_b = A_m(entry1,1);
B_b = A_m(entry1,2);
C_b = A_m(entry2,1);
D_b = A_m(entry2,2);
Z_b(j) = (A_b*Z_NB + B_b)/(C_b*Z_NB + D_b);
end

%% 2x2 matrices for parallel connection of piezoelectric layers

A_c = zeros(length(omega),1);
B_c = zeros(length(omega),1);
C_c = zeros(length(omega),1);
D_c = zeros(length(omega),1);

for jj = 1:length(omega)
A_c(jj) = T_31 - T_33*(T_21(jj)*Z_b(jj) + T_11(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
B_c(jj) = T_32 - T_33*(T_22(jj)*Z_b(jj) + T_12(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
C_c(jj) = T_41(jj) - T_43(jj)*(T_21(jj)*Z_b(jj) + T_11(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
D_c(jj) = T_42(jj) - T_43(jj)*(T_22(jj)*Z_b(jj) + T_12(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
end

%A_C = zeros(2*length(omega),2);

for x = 1:length(omega)
    entry1 = 2*x - 1;
    entry2 = 2*x;
A_C(entry1:entry2,1:2) = [A_c(x) B_c(x); C_c(x) D_c(x)]^2;
end

%% Final 2x2 matrix of Piezo stack

A_ps = A_C;

A = zeros(length(omega),1);
B = zeros(length(omega),1);
C = zeros(length(omega),1);
D = zeros(length(omega),1);



for xx = 1:length(omega)
    entry1 = 2*xx - 1;
    entry2 = 2*xx;
AA_ps = A_C(entry1,1);
BB_ps = A_C(entry1,2);
CC_ps = A_C(entry2,1);
DD_ps = A_C(entry2,2);

A_ps = [AA_ps BB_ps; CC_ps DD_ps];

A(xx) = A_ps(1,1);
B(xx) = A_ps(1,2);
C(xx) = A_ps(2,1);
D(xx) = A_ps(2,2);
end


%% Sensitivity Functions

Z_E = zeros(length(omega),1);
V_IL = zeros(length(omega),1);

for a = 1:length(omega)
Z_E(a) = (A(a)*Z_FP + B(a))/(C(a)*Z_FP + D(a));            % electrical input impedance
V_IL(a) = Z_FP/(A(a)*Z_FP + B(a));                         % voltage transfer ratio between input voltage and output force
end

phase_impedance = rad2deg(angle(Z_E));
phase_tf = rad2deg(angle(V_IL));
mag_impedance = 20*log10(abs(Z_E));
mag_tf = 20*log10(abs(V_IL));

[pks1 ind1] = findpeaks(mag_impedance,'MinPeakDistance',5000,'MinPeakProminence',0.01);
figure
subplot(2,1,1); plot(f(ind1),pks1,'or'); hold on; plot(f,mag_impedance); title('Bode magnitude plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');
subplot(2,1,2); plot(f,phase_impedance); title('Phase plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');


[pks ind] = findpeaks(mag_tf,'MinPeakDistance',5000,'MinPeakProminence',0.01);

figure
subplot(2,1,1); plot(f(ind),pks,'or'); hold on; plot(f,mag_tf); title('Bode magnitude plot of the transfer function V_{in} --> F_{out}'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');
subplot(2,1,2); plot(f,phase_tf); title('Phase plot of the transfer function V_{in} --> F_{out}'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');






    
    