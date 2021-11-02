function [v] = EOM_Particle(Fac)
%EOM_PARTICLE Summary of this function goes here
%   Detailed explanation goes here
syms seis pyrs r mu
v = Fac/(seis*pyrs*r*mu);
end

