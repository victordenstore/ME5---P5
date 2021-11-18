close all;
clear all;
clc;

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
Z_in_mag = 20*log10(abs(Z_in))

figure
tiledlayout(2,1); nexttile
plot(M(:,1),Z_in_mag); title('Impedance magnitude'); xlabel('frequency [kHz]'); ylabel('magnitude [DB]');
nexttile; plot(M(:,1),-phase1); title('Impedance angle'); xlabel('frequency [kHz]'); ylabel('angle [degree]');

 