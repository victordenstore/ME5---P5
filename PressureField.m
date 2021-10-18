function P_local = PressureField(w,t,x,lambda,damp_coeff,P)
%PRESSUREFIELD Summary of this function goes here
%   Detailed explanation goes here


    P_local = P * exp(1i*w*t-1i*x*2*pi/lambda - damp_coeff * x);

end

