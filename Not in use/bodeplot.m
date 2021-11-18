close all;
clear all;
clc;

x=0.0000042;

%% intialization
f = 1*10^6;
w=linspace(1,4*pi*f,10000); % System frequency in rad/s, w=linspace(2.3064*pi*f,2.3068*pi*f,10000); 
a=2.5*10^-3;                     % Radius of the piezoelectric plate
S=pi*a^2;                        % Piezoelectric surface area
rho=7.8*10^3;                    % Piezoelectric plate density
c_33=16.6*10^10;                 % Elastic constant of the plate
d=2*10^-3;                       % plate thickness
rho_4=7860;                          % Density of plain carbon steel 
c_3=3230;                            % Speed of wave in plain carbon steel
Z_b=rho_4*c_3*S;                         % Acoustic impedance of the backing plate
d_33=265*10^-12;                 % Piezoelectric charge coefficient
s_33=14.2*10^-12;                % Elastic compliance coefficient
epsilon_0=8.854*10^-12;          % Vacuum permitivity
epsilon_33=1200*epsilon_0;       % relativity Permativity
Beta_33= epsilon_0/epsilon_33;   % the dielectric impermeability of the plate at constant strain,
c=1300;                          % compresisonal wave speed in fluid                     
h_33= d_33/(s_33*epsilon_33);    % Piezoelectric stiffness constant for the plate
v_o=sqrt(c_33/rho);              % compressional wave speed from piezoelectric
k=w/v_o;                         % wave number for the peizoelectric plate
C_o=2*S*epsilon_33/(2*10^-3);               % the clamped capacitance of the plate
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
d_2=0.5*10^-3;                    % thickness of the coating material of the transducer source
S_m=(4.5*10^-3)^2*pi;             % Face area of the coating material fo the transducer soruce
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
multiply1(i) = 1 / (Z_b - 1i*Z_o*tan(k(i)*d/2));
T_A_11(i) = (Z_b + 1i*Z_o*cot(k(i)*d))/n + n/(1i*w(i)*C_o);
T_A_12(i) = ((Z_o)^2+1i*Z_o*Z_b*cot(k(i)*d))/n + n*(Z_b - 2*1i*Z_o*tan(k(i)*d/2))/(1i*w(i)*C_o); 
T_A_21(i) = -1i*w(i)*C_o*(Z_b+1i*Z_o*cot(k(i)*d));
T_A_22(i) = -1i*w(i)*C_o*((Z_o)^2 + 1i*Z_o*Z_b*cot(k(i)*d));
T_A(boy:boy1,1:2) = multiply1(i) * [T_A_11(i) T_A_12(i); T_A_21(i) T_A_22(i)];
TAl_matrix(boy:boy1,1:2)=[cos(k_2(i)*d_2) -1i*Z_m*sin(k_2(i)*d_2); (-1i*sin(k_2(i)*d_2))/Z_m cos(k_2(i)*d_2)]; 

Z_Aa_r = rho_2*c*S_a;
T_A = T_A(boy:boy1,1:2)*TAl_matrix(boy:boy1,1:2);
S_A_vl(i) = 1/(Z_Aa_r*T_A(2,1) + T_A(2,2));

Z_Ae_in(i) = (Z_Aa_r*T_A(1,1) + T_A(1,2))/(Z_Aa_r*T_A(2,1) + T_A(2,2))+1i*w(i)*4.34*10^-6;

S_FV(i) = Z_Aa_r * S_A_vl(i)/Z_Ae_in(i);
end

phase = rad2deg(angle(S_FV));

S_FV = 20*log10(abs(S_FV))
%% figure

figure
tiledlayout(2,1);
nexttile
plot(w/(2*pi),S_FV,'g')
ylim([-75 -25])
nexttile
plot(w/(2*pi),phase)

figure
tiledlayout(2,1);
nexttile
Z_Ae_in1 =  20*log10(abs(Z_Ae_in))
plot(w/(2*pi),Z_Ae_in1,'g')
nexttile
phase1 = rad2deg(angle(Z_Ae_in));
plot(w/(2*pi),phase1)


M = readmatrix('C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\FreqSweepSkabelon.xlsx')

V_in = 75;
I_phasor = zeros(length(M),1); 
Z_in = zeros(length(M),1);
      
for i = 1:length(M)
    phase = exp(1i*M(i,4)*pi/180);
    I_phasor(i) = M(i,3)*phase;
end
      for j = 1:length(M)
          Z_in(j) = V_in/I_phasor(j);
      end

phase1 = rad2deg(angle(Z_in));
Z_in_mag = 20*log10(abs(Z_in));

figure
tiledlayout(2,1); nexttile
plot(M(:,1),Z_in_mag); title('Impedance magnitude'); xlabel('frequency [kHz]'); ylabel('magnitude [DB]');
nexttile; plot(M(:,1),-phase1); title('Impedance angle'); xlabel('frequency [kHz]'); ylabel('angle [degree]');


%% 2nd transducer, pressure wave magnitude calculation

R1 = 4.2736;                                                                % Left wall reflection coefficient.
R2 = 4.8295;                                                                % Right wall reflection coefficient.

w = linspace(100e3*2*pi,1500e3*2*pi,141);                                               % time duration during which the pressure field
                                                                            % is monitored.
x = 4.48e-2;                                                                % Distance over which the pressure field is monitored.
P_field = zeros(length(w),1);
c = 1300;

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
multiply1(i) = 1 / (Z_b - 1i*Z_o*tan(k(i)*d/2));
T_A_11(i) = (Z_b + 1i*Z_o*cot(k(i)*d))/n + n/(1i*w(i)*C_o);
T_A_12(i) = ((Z_o)^2+1i*Z_o*Z_b*cot(k(i)*d))/n + n*(Z_b - 2*1i*Z_o*tan(k(i)*d/2))/(1i*w(i)*C_o); 
T_A_21(i) = -1i*w(i)*C_o*(Z_b+1i*Z_o*cot(k(i)*d));
T_A_22(i) = -1i*w(i)*C_o*((Z_o)^2 + 1i*Z_o*Z_b*cot(k(i)*d));
T_A(boy:boy1,1:2) = multiply1(i) * [T_A_11(i) T_A_12(i); T_A_21(i) T_A_22(i)];
TAl_matrix(boy:boy1,1:2)=[cos(k_2(i)*d_2) -1i*Z_m*sin(k_2(i)*d_2); (-1i*sin(k_2(i)*d_2))/Z_m cos(k_2(i)*d_2)]; 

Z_Aa_r = rho_2*c*S_a;
T_A = T_A(boy:boy1,1:2)*TAl_matrix(boy:boy1,1:2);
S_A_vl(i) = 1/(Z_Aa_r*T_A(2,1) + T_A(2,2));

Z_Ae_in(i) = (Z_Aa_r*T_A(1,1) + T_A(1,2))/(Z_Aa_r*T_A(2,1) + T_A(2,2))+1i*w(i)*4.34*10^-6;

S_FV(i) = Z_Aa_r * S_A_vl(i)/Z_Ae_in(i);
end

V_M = 75;                                                                    % Input voltage
F_surface = V_M * S_FV;                                                       % Total force on trasnducer surface
P_surface = F_surface/(S_a);

for n = 1:length(w)
    w1 = w(n);
    f = w1/(2*pi);
    T = 1/f;
    t = T/4;
    lambda = c/f;
     
     P_field(n) = P_surface(n)*(R1+1) * exp(1i*w1*t-1i*x*2*pi/lambda) ...
         - R2*P_surface(n) * exp(1i*w1*t+1i*x*2*pi/lambda);
      end
 
F_field = P_field*S;
volt_trans_2 = F_field./S_FV;
mag_PF = 20*log10(abs(volt_trans_2));
freq = w/(2000*pi);

figure
subplot(2,1,1)
plot(freq,mag_PF); title('Theoretical voltage at 2nd transducer'); xlabel('frequency [kHz]'); ylabel('magnitude [V]');
subplot(2,1,2); plot(M(:,1),M(:,5)); title('2nd transducer magnitude'); xlabel('frequency [kHz]'); ylabel('magnitude [mV]');

%% Fourier domain, of the sine sweep

filename = 'C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\FrequencySweep400k1.5M40ms_10ms_000_ALL.csv';
sinesweep = readtable(filename);
n = 617385:1255020;
%n = 1:length(sinesweep.CH2);
Fs = 1/2e-8;                    % Note that the sampling time actually fluctuates, but idk how to handle that
T = 1/Fs;
L = length(n);
t = (0:L-1)*T;


current = sinesweep.CH2(n);
voltage = sinesweep.CH3(n);
impedance = voltage./current;

inf_values = find(isinf(impedance));
newvalues = zeros(length(inf_values),1);

for a = 1:length(inf_values)
    newvalues(a) = (impedance(inf_values(a)+1) - impedance(inf_values(a)-1))/2;
end

impedance(inf_values) = newvalues;


timelength = sinesweep.TIME(n);
figure
subplot(3,1,1); plot(timelength,current) 
title('current in the time domain')
xlabel('time [s]')
ylabel('Amplitude [v]');
subplot(3,1,2); plot(timelength,voltage) 
title('voltage in the time domain')
xlabel('time [s]')
ylabel('Amplitude [v]');
subplot(3,1,3); plot(timelength,impedance) 
title('impedance in the time domain')
xlabel('time [s]')
ylabel('Amplitude [v]');

figure
plot(timelength,impedance)


lowpass(impedance,1.5e6,Fs)
Y = fft(current);
Y1 = fft(voltage);
Y2 = fft(impedance);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P4 = abs(Y1/L);
P3 = P4(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);

P6 = abs(Y2/L);
P5 = P6(1:L/2+1);
P5(2:end-1) = 2*P5(2:end-1);


f = Fs*(0:(L/2))/L;
figure
subplot(3,1,1); plot(f,P1) 
title('Single-Sided Amplitude Spectrum of the current')
xlabel('f (Hz)')
subplot(3,1,2); plot(f,P3) 
title('Single-Sided Amplitude Spectrum of the voltage')
xlabel('f (Hz)')
ylabel('|P1(f)|')
subplot(3,1,3); plot(f,P5) 
title('Single-Sided Amplitude Spectrum of the impedance')
xlabel('f (Hz)')
ylabel('|P1(f)|')




