close all;
clear all;
clc;

close all;
clear all;
clc;

f = [1 2 3 4 5 6 7 8 8.3 8.6 8.8 9 9.2 9.5 9.8 10 10.2 11 12 13 14 15]*10^5;
V_in = 75;
I_in = [11.5 23.19 32.62 49.21 40.4 54.42 65.22 74.41 65.9228 65.8047 65.7044 79.15 ...
         65.3308 64.9142 64.4894 91.93 64.0824 95.25 92.81 98.35 107.3 119.9];
I_in = I_in/2 * 10^-3;
theta = [78.72 68.44 63.94 50.93 41.44 42.08 35.34 26.87 ...
         24.9276 23.2507 22.2266 21.24  20.2422 18.4695 15.9735 13.69 10.8722 4.16 -1.1 -1.4 348.1 348]';
I_phasor = zeros(length(f),1); 
Z_in = zeros(length(f),1);


      
for i = 1:length(f)
    phase = exp(1i*theta(i)*pi/180);
    I_phasor(i) = I_in(i)*phase;
end
      for j = 1:length(f)
          Z_in(j) = V_in/I_phasor(j);
      end

phase1 = rad2deg(angle(Z_in));
Z_in_mag = 20*log10(abs(Z_in))

figure
tiledlayout(2,1); nexttile
plot(f,Z_in_mag)
nexttile; plot(f,phase1)

 