close all; clear all; clc;


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

phase_impedance = rad2deg(angle(impedance_phasor));
mag_impedance = 20*log10(abs(impedance_phasor));

admittance = 1./impedance_phasor;
phase_admittance = rad2deg(angle(admittance));
mag_admittance = 20*log10(abs(admittance));
susceptance = imag(admittance);

figure
subplot(2,1,1); plot(freq_vec,mag_impedance); title('Bode magnitude plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]'); ylim([40 90])
subplot(2,1,2); plot(freq_vec,phase_impedance); title('Phase plot of the impedance'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');

ptp_trans2_vec = ptp_trans2_vec./2000;
mag_trans2 = 20*log10(abs(ptp_trans2_vec));

figure
subplot(2,1,1); plot(freq_vec,ptp_trans2_vec); title('Bode magnitude plot of the second transducer'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [V]');

figure
subplot(2,1,1); plot(freq_vec,mag_admittance); title('Bode magnitude plot of the admittance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');
subplot(2,1,2); plot(freq_vec,phase_admittance); title('Phase plot of the admittance'); ...
    xlabel('frequency [Hz]'); ylabel('phase [degree]');

figure
plot(freq_vec,susceptance); title('Bode magnitude plot of the susceptance'); ...
    xlabel('frequency [Hz]'); ylabel('magnitude [DB]');


