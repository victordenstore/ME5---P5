close all;
clear all;
clc;

x=0.0000042;

%% intialization
f = 1*10^6;
w=linspace(1,2.5*pi*f,100000);                     % System frequency in rad/s    
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
L=0.00448;   

c_2 = 3230;                      % wave speed in steel
k_2=w/c_2                       % wave number of the coating material of the transducer source
d_2=0.5*10^3;                    % thickness of the coating material of the transducer source
S_m=(4.5*10^3)^2*pi;             % Face area of the coating material fo the transducer soruce
Z_m=rho_3*c_2*S_m;               % The acoustic impedance of the coating material of the transducer source
%% TA Matrix
T_A_11 = zeros(length(w),1);
T_A_12 = zeros(length(w),1);
T_A_21 = zeros(length(w),1);
T_A_22 = zeros(length(w),1);
multiply1 = zeros(length(w),1);
T_A = zeros(2*length(w),2);
S_A_vl = zeros(length(w),1);
Z_Ae_in = zeros(length(w),1);
S_FV  = zeros(length(w),1);




for i = 1:length(w)
    
    boy = 2*i-1;
    boy1 = 2*i;
multiply1(i) = 1 / (Z_b(i) - 1i*Z_o*tan(k(i)*d/2));
T_A_11(i) = (Z_b(i) + 1i*Z_o*cot(k(i)*d))/n + n/(1i*w(i)*C_o);
T_A_12(i) = ((Z_o)^2+1i*Z_o*Z_b(i)*cot(k(i)*d))/n + n*(Z_b(i) - 2*1i*Z_o*tan(k(i)*d/2))/(1i*w(i)*C_o); 
T_A_21(i) = -1i*w(i)*C_o*(Z_b(i)+1i*Z_o*cot(k(i)*d));
T_A_22(i) = -1i*w(i)*C_o*((Z_o)^2 + 1i*Z_o*Z_b(i)*cot(k(i)*d));
T_A(boy:boy1,1:2) = multiply1(i) * [T_A_11(i) T_A_12(i); T_A_21(i) T_A_22(i)];
TAl_matrix(boy:boy1,1:2)=[cos(k_2(i)*d_2) -1i*Z_m*sin(k_2(i)*d_2); (-1i*sin(k_2(i)*d_2))/Z_m cos(k_2(i)*d_2)]; 

Z_Aa_r = rho_2*c*S_a;
T_A = T_A(boy:boy1,1:2)*TAl_matrix(boy:boy1,1:2);
S_A_vl(i) = 1/(Z_Aa_r*T_A(2,1) + T_A(2,2));

Z_Ae_in(i) = (Z_Aa_r*T_A(1,1) + T_A(1,2))/(Z_Aa_r*T_A(2,1) + T_A(2,2));

S_FV(i) = Z_Aa_r * S_A_vl(i)/Z_Ae_in(i);
end

S_FV = 20*log10(abs(S_FV))
%% figure

tiledlayout(2,1);
nexttile
plot(w,S_FV)
phase = rad2deg(angle(S_FV));
nexttile
plot(w,phase)

