x=0.0000042;

% variables
w=1*10^6;                  % System frequency in rad/s
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

V_in=150*exp(1i*w);     %The input voltage (frequency domain)

F=V_in*S_FV;            %The force output with a certain frequency
P=abs(F)/S;             %Pressure at the transducer face