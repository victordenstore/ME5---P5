% function Fac = Gorkov(acoustic_contrast_factor,r,Eac,k,x)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% Fac = 4*pi*acoustic_contrast_factor*r*3*k*Eac*sin(k*x);
% end
         function Fac = Gorkov(c, rho_l, r, f_1, f_2,w,T,x,damp_coeff,P,k)


grad_p_in = P*((1i*w*T - 1)*exp(1i*w*T) + 1)/(1i*w*T^2)...
     - damp_coeff * P*exp(-damp_coeff*x)/(1i*w*T) * (exp(1i*w*T)-1);
% grad_p_in = P*exp(-damp_coeff*x)*(exp(1i*w*T)-1)/(1i*w*T);
grad_v_in = P*((1i*w*T - 1)*exp(1i*w*T) + 1)/(1i*w*T^2*rho_l*c)...
    - damp_coeff * P*exp(-damp_coeff*x)/(1i*w*T*rho_l*c) * (exp(1i*w*T)-1);

Fac = 4*pi/3 * r^3 * (f_1 * 1/(2*rho_l*c^2) * grad_p_in - f_2 * 3/4*rho_l*grad_v_in);


end