clear all;
close all;
clc;

filename = 'C:\Users\wneum\OneDrive - Universiteit Twente\Desktop\project\FrequencySweep1_000_ALL.csv';
M = readtable(filename);




figure
subplot(2,1,1), plot(M.TIME,M.CH2); 
subplot(2,1,2), plot(M.TIME,M.CH3);




