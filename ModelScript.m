clear all; close all; clc;
x=0.0000042;

% variables
w=2*pi*10^6;                  % System frequency in rad/s
a=2.5*10^-3;                     % Radius of the piezoelectric plate
rho=7.8*10^3;                    % Piezoelectric plate density
c_33=16.6*10^10;                 % Elastic constant of the plate
d=2*10^-3;                       % plate thickness
Z_b=x*w;                        % Acoustic impedance of the backing plate (this is a function of frequency)
d_33=265*10^-12;                 % Piezoelectric charge coefficient
s_33=14.2*10^-12;                % Elastic compliance coefficient
epsilon_0=8.854*10^-12;          % Vacuum permitivity
epsilon_33=1200*epsilon_0;       % relativity Permativity
Beta_33= epsilon_0/epsilon_33;   % the dielectric impermeability of the plate at constant strain,
% equations
h_33= d_33/(s_33*epsilon_33);    % Piezoelectric stiffness constant for the plate
v_o=sqrt(c_33/rho);              % compressional wave speed from piezoelectric
k=w/v_o;                         % wave number for the peizoelectric plate
S=pi*a^2;                        % Piezoelectric surface area
C_o=S/(Beta_33*d);               % the clamped capacitance of the plate
n=h_33*C_o;                      % A given constant
Z_o=rho*v_o*S;                   % plane wave acoustic impedance of the piezoelectric plate

% Matricies
TAe_matrix=[1/n n/(1i*w*C_o); -1i*w*C_o 0];          %Transducer electrical matrix
TAa_matrix=(1/(Z_b-1i*Z_o*tan(k*d/2)))*[Z_b+1i*Z_o*cot(k*d) (Z_o)^2+1i*Z_o*Z_b*cot(k*d); 1 Z_b-2*1i*Z_o*tan(k*d/2)];     %Transducer acoustic matrix

TA_matrix=TAe_matrix*TAa_matrix; %Sittig model matrix

%equations
c=1500;                  % compresisonal wave speed in fluid
S_a=S;                   % effective face area of the transducer
rho_2=900;               % density of the fluid
Z_r=S_a*c*rho_2;         % Acoustic radiation impedance - just in case we want to put in the Ka condition later.

S_vI=1/(Z_r*TA_matrix(2,1)+TA_matrix(2,2)); %Sensitivity for vI.

Z_in=(Z_r*TA_matrix(1,1)+TA_matrix(1,2))/(Z_r*TA_matrix(2,1)+TA_matrix(2,2)); %Electrical impedance of the transducer

S_FV=Z_r*S_vI/Z_in; %Sensitivity FV


t = linspace(0,10*pi/w,1000);
x = linspace(0,50*10^-3,1000);
alfa = -40;
V_in = zeros(1000,1);
F = zeros(1000,1);
P = zeros(1000,1);
expression = zeros(1000,1000);

% for ko = 1:1000
%     for lim = 1:1000
% V_in(ko) = 150*exp(1i*w*t(ko));     %The input voltage (frequency domain/time)
% V_in(ko) = 150*sin(w*t(ko));     %The input voltage (frequency domain/time)
% F(ko) =  V_in(ko)*S_FV;            %The force output with a certain frequency
% P(ko) = F(ko)/S;                   %Pressure at the transducer face
% expression(ko,lim) = (P(ko)) * exp(alfa*x(lim));
%     end
% end
V_M = 150;
lambda = 1.5*10^-3;
P = S_FV * V_M/S;
damp_coeff = 60;
exp1 = zeros(1000,1000);

for to = 1:1000
    for gi = 1:1000
    exp1(to,gi) = P * exp(1i*w*t(gi)-1i*x(to)*2*pi/lambda - damp_coeff * x(to));
    end
end







% [x,t] = meshgrid(x,t);
% surf(x,t,abs(expression)

figure
contourf(t, x, real(exp1), 14); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); 
figure
meshc(t, x, real(exp1)); colormap jet; colorbar;
xlabel('time [s]'); ylabel('distance [m]'); zlabel('pressure [Pa]')

s = tf('s');     %initialize laplace variable
rho_air = 1.225; %kg/m^3
rho_oil = 0.91;
visc_oil = 1.3; 
R = 0.05*10^-3;
A_b = pi*R^2;
m_b = rho_air * (4/3) * pi * R^3;
tf_particle = A_b/(m_b * (1i*w)^2 + 6*pi*R*rho_oil*visc_oil*(1i*w));

p_grad = zeros(1000,1000);

for ji = 2:999
    for hi = 1:1000
        p_grad(hi,ji) = exp1(hi,ji+1) - exp1(hi,ji-1);
    end
end

tf_Vin_xout = tf_particle * p_grad;

x_0 = 300;
t_0 = 100;
x1 = zeros(300,1);
t = linspace(0,10*pi/w,1000);

time_vec = linspace(0,10*pi/(w*(1000/3.3)),300);

for gif = 1:300
    for dsf = 1:300
    x1(gif) = p_grad(t_0+dsf,x_0+gif) * tf_particle;
    end
end

figure
for hello = 1:300
plot(time_vec(hello), real(x1(hello)))
hold on
plot(time_vec(1:hello), real(x1(1:hello)))

xlim([0 2*10^-8])
ylim([-1*10^-8 1*10^-8])
pause(0.05);
end
