close all;
clear all;
clc;
%% Parameter initialization

j = 1i;                                                                     % imaginary value
f = linspace(10000,2*10^6,10000);                                                               % operating frequency [Hz]
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


gamma = pi*f*d_n;                                                    % 
eps_0=8.854*10^-12;          
eps_33=1200*eps_0;
s_33 = 14.2*10^-12;
d_33 = 265*10^-12;
h_33 = d_33/(s_33*eps_33);

k_n_squared = 0.66^2;
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


% T_11 = (cos(gamma) - s)./(1-s);
% T_12 = (j*Z_0*(sin(gamma) - 2*c))./(1-s);
% cosinus_gamma = cos(gamma)-1;
% T_13 = -cosinus_gamma.*phi./(1-s);
% T_14 = 0;
% T_21 = (j*sin(gamma))./(Z_0*(1-s));
% T_22 = (cos(gamma) - s)./(1-s);
% T_23 = -(j*phi.*sin(gamma))./(1-s);
% T_24 = 0;
% T_31 = 0;
% T_32 = 0;
% T_33 = 1;
% T_34 = 0;
% T_41_intermediate = -(j*sin(gamma))./(Z_0*(1-s));
% T_41 = T_41_intermediate .* phi;
% T_42 = -((cos(gamma) - 1).*phi)./(1-s);
% T_43 = (j*omega*C_0)./(1-s);
% T_44 = 1;
% 
% Tn = zeros(4*length(f),4);
% 
% for hh = 1:length(f)
%     row1 = 4*hh - 3;
%     row2 = 4*hh - 2;
%     row3 = 4*hh - 1;
%     row4 = 4*hh;
%     rows = [row1 row2 row3 row4];
%     Tn(rows(1):rows(4),1:4) =[T_11(hh) T_12(hh) T_13(hh) T_14; ...
%                     T_21(hh) T_22(hh) T_23(hh) T_24; ...
%                     T_31 T_32 T_33 T_34; ...
%                     T_41(hh) T_42(hh) T_43(hh) T_44];
% end
% 
% r_copper = 2.5*10^-3;
% A_copper = pi*r_copper^2;
% rho_copper = 8920;
% c_copper = 3570;
% d_copper = 0.5*10^-3;
% v_copper = sqrt(c_copper/rho_copper);
% 
% Z_copper = A_copper*rho_copper*c_copper;
% gamma_copper = 0;
% 
% A_copper = [cos(gamma_copper) j*Z_copper*sin(gamma_copper); ...
%     (j*sin(gamma_copper))/Z_copper cos(gamma_copper)];
% A_E = 1;
% B_E = 0;
% C_E = 0;
% D_E = 1;
% 
% T_intermediate_copper = [A_copper(1,1) A_copper(1,2) 0 0; ...
%                   A_copper(2,1) A_copper(2,2) 0 0; ...
%                   0 0 A_E B_E; ...
%                   0 0 C_E D_E];
%                             
% rho_FP = 7850;
% r_FP = 4.5*10^-3;
% c_FP = 3230;
% S_FP = pi*r_FP^2;
% d_FP = 0.5*10^-3;
% v_FP = sqrt(c_FP/rho_FP);
% 
% Z_FP = rho_FP*c_FP*S_FP;
% gamma_FP = pi*f*d_FP;
% 
% A_E_FP = 1;
% B_E_FP = 1./(omega*1400*eps_0*S_FP/d_FP);
% C_E_FP = 0;
% D_E_FP = 1;
% 
% A_FP = zeros(2*length(f),2);
% T_intermediate_FP = zeros(4*length(f),4);
% 
% for ee = 1:length(f)
%     row1 = 4*ee - 3;
%     row2 = 4*ee - 2;
%     row3 = 4*ee - 1;
%     row4 = 4*ee;
%     rows = [row1 row2 row3 row4];
% A_FP(rows(1):rows(2),1:2) = [cos(gamma_FP(ee)) j*Z_FP*sin(gamma_FP(ee)); ...
%     (j*sin(gamma_FP(ee)))/Z_FP cos(gamma_FP(ee))];
% A_FP_now = A_FP(rows(1):rows(2),1:2);
% T_intermediate_FP(rows(1):rows(4),1:4) = [A_FP_now(1,1) A_FP_now(1,2) 0 0; ...
%                     A_FP_now(2,1) A_FP_now(2,2) 0 0; ...
%                     0 0 A_E_FP B_E_FP(ee);
%                     0 0 C_E_FP D_E_FP];
% end
% 
% T = zeros(4*length(f),4);
% T_11 = zeros(length(f),1);
% T_12 = zeros(length(f),1);
% T_13 = zeros(length(f),1);
% T_14 = zeros(length(f),1);
% T_21 = zeros(length(f),1);
% T_22 = zeros(length(f),1);
% T_23 = zeros(length(f),1);
% T_24 = zeros(length(f),1);
% T_31 = zeros(length(f),1);
% T_32 = zeros(length(f),1);
% T_33 = zeros(length(f),1);
% T_34 = zeros(length(f),1);
% T_41 = zeros(length(f),1);
% T_42 = zeros(length(f),1);
% T_43 = zeros(length(f),1);
% T_44 = zeros(length(f),1);
% 
% for tt = 1:length(f)
%     row1 = 4*tt - 3;
%     row2 = 4*tt - 2;
%     row3 = 4*tt - 1;
%     row4 = 4*tt;
%     rows = [row1 row2 row3 row4];
% T(rows(1):rows(4),1:4) = Tn(rows(1):rows(4),1:4)*Tn(rows(1):rows(4),1:4)* T_intermediate_copper^2 * T_intermediate_FP(rows(1):rows(4),1:4);
% T = T(rows(1):rows(4),1:4);
% T_11(tt)=Tn(1,1);
% T_12(tt)=Tn(1,2);
% T_13(tt)=Tn(1,3);
% T_14(tt)=Tn(1,4);
% T_21(tt)=Tn(2,1);
% T_22(tt)=Tn(2,2);
% T_23(tt)=Tn(2,3);
% T_24(tt)=Tn(2,4);
% T_31(tt)=Tn(3,1);
% T_32(tt)=Tn(3,2);
% T_33(tt)=Tn(3,3);
% T_34(tt)=Tn(3,4);
% T_41(tt)=Tn(4,1);
% T_42(tt)=Tn(4,2);
% T_43(tt)=Tn(4,3);
% T_44(tt)=Tn(4,4);
% 
% end

%% parameters

f = linspace(1,2*10^6,10000);
omega = 2*pi*f;
N = 2;
j = 1i;
C_s = S*eps_33/d_n;
h = h_33;
Z = S*rho_n*v_n;
Z_inv = 1/Z;
theta = omega*d_n/v_n;
k = 0.66;
sigma = k^2./theta;
phi = acos((cos(theta)-sigma.*sin(theta))./(1-sigma.*sin(theta)));
R = (sqrt(sin(theta)-2*sigma.*(1-cos(theta))))./sin(theta);
R_inv = 1./R;



T_11 = cos(N*phi);
T_12 = -j*Z*R.*sin(N*phi);
T_13 = -h*C_s*tan(0.5*phi).*sin(N*phi);
T_14 = 0;
T_21 = -j*Z_inv*R_inv.*sin(N*phi);
T_22 = cos(N*phi);
T_23 = -j*h*C_s*Z_inv*R_inv.*tan(0.5*phi).*(cos(N*phi)-(-1)^N);
T_24 = 0;
T_31 = 0;
T_32 = 0;
T_33 = (-1)^N;
T_34 = 0;
T_41 = -j*h*C_s*Z_inv*R_inv.*tan(0.5*phi).*(cos(N*phi)-(-1)^N);
T_42 = -h*C_s*tan(0.5*phi).*sin(N*phi);
T_43 = j*((N*(-1)^N)*(1+2*sigma.*R_inv.*tan(0.5*phi))+sigma.*R_inv.*((tan(0.5*phi)).^2).*sin(N*phi)).*omega*C_s;
T_44 = (-1)^N;


%% Front and back layers (loaded impedance)

rho_H = 7860;
r_H = 2.5*10^-3;
c_H = 3230;
S_H = pi*r_H^2;
d_H = 3.25*10^-3;
v_H = sqrt(c_H/rho_H);
gamma_H = pi*f*d_H;
Z_H = rho_H*c_H*S_H;

A_H = zeros(2*length(f),2);

for i = 1:length(omega)
    entry1 = 2*i - 1;
    entry2 = 2*i;
A_H(entry1:entry2,1:2) = [cos(gamma_H(i)) j*Z_H*sin(gamma_H(i)); ...
    (j*sin(gamma_H(i)))/Z_H cos(gamma_H(i))];
end

rho_B = 7860;
r_B = 2.5*10^-3;
c_B = 3230;
S_B = pi*r_B^2;
d_B = 20*10^-3;
v_B = sqrt(c_H/rho_H);
gamma_B = pi*f*d_B;
Z_B = rho_H*c_H*S_H;

A_B = zeros(2*length(f),2);

for ii = 1:length(omega)
    entry1 = 2*ii - 1;
    entry2 = 2*ii;
A_B(entry1:entry2,1:2) = [cos(gamma_B(ii)) j*Z_B*sin(gamma_B(ii)); ...
    (j*sin(gamma_B(ii)))/Z_B cos(gamma_B(ii))];
end

Z_b = zeros(length(f),1);
A_m = zeros(2*length(f),2);

for j = 1:length(f)
    entry1 = 2*j - 1;
    entry2 = 2*j;
    entries = [entry1 entry2];
A_H1 = A_H(entries(1):entries(2),1:2);
A_B1 = A_B(entries(1):entries(2),1:2);
A_m(entries(1):entries(2),1:2) = A_H1*A_B1;
A_b = A_B(entry1,1);
B_b = A_B(entry1,2);
C_b = A_B(entry2,1);
D_b = A_B(entry2,2);
Z_b(j) = (A_b*Z_B + B_b)/(C_b*Z_B + D_b);
end

%% 2x2 matrices for parallel connection of piezoelectric layers

A_c = zeros(length(f),1);
B_c = zeros(length(f),1);
C_c = zeros(length(f),1);
D_c = zeros(length(f),1);

for jj = 1:length(omega)
A_c(jj) = T_31 - T_33*(T_21(jj)*Z_b(jj) + T_11(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
B_c(jj) = T_32 - T_33*(T_22(jj)*Z_b(jj) + T_12(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
C_c(jj) = T_41(jj) - T_43(jj)*(T_21(jj)*Z_b(jj) + T_11(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
D_c(jj) = T_42(jj) - T_43(jj)*(T_22(jj)*Z_b(jj) + T_12(jj))/(T_23(jj)*Z_b(jj) + T_13(jj));
end

A_C = zeros(2*length(f),2);

for x = 1:length(f)
    entry1 = 2*x - 1;
    entry2 = 2*x;
A_C(entry1:entry2,1:2) = [A_c(x) B_c(x); C_c(x) D_c(x)];
end

%% Final 2x2 matrix of Piezo stack

A_ps = zeros(2*length(f),2);
A = zeros(length(f),1);
B = zeros(length(f),1);
C = zeros(length(f),1);
D = zeros(length(f),1);

for rr = 1:length(f)
entry1 = 2*rr - 1;
entry2 = 2*rr;

A_ps(entry1:entry2,1:2) = A_C(entry1:entry2,1:2)*A_m(entry1:entry2,1:2);
A_ps_now = A_ps(entry1:entry2,1:2);
A(rr) = A_ps_now(1,1);
B(rr) = A_ps_now(1,2);
C(rr) = A_ps_now(2,1);
D(rr) = A_ps_now(2,2);
end


%% Sensitivity Functions

Z_E = zeros(length(omega),1);
V_IL = zeros(length(omega),1);

for a = 1:length(f)
Z_E(a) = (A(a)*Z_H + B(a))/(C(a)*Z_H + D(a));            % electrical input impedance
V_IL(a) = Z_H/(A(a)*Z_H + B(a));                         % voltage transfer ratio between input voltage and output force
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







    
    