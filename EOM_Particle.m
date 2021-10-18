function [a] = EOM_Particle(Mp,Fac,B,v)
%EOM_PARTICLE Summary of this function goes here
%   Detailed explanation goes here
a = 1/Mp * (Fac - B*v);
end

