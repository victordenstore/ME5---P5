clear all;
close all;
clc;

filename = 'C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\FrequencySweep400k1.5M40ms_10ms_000_ALL.csv';
M = readtable(filename);

n = 617385:1255020;

figure
subplot(2,1,1); plot(M.TIME(n),M.CH2(n));
subplot(2,1,2); plot(M.TIME(n),M.CH3(n));


impedance = M.CH2(n)./M.CH3(n);
inf_values = find(isinf(impedance));
newvalues = zeros(length(inf_values),1);

for a = 1:length(inf_values)
    newvalues(a) = (impedance(inf_values(a)+1) - impedance(inf_values(a)-1))/2;
end

impedance(inf_values) = newvalues;

figure
plot(M.TIME(n),impedance);

g = spa(M.CH2);
p = etfe(M.CH2);
figure
spectrum(g,p);