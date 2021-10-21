         function Fac = Gorkov(c, rho_l, r, f_1, f_2,w,T,x,damp_coeff,P)


grad_p_in = P*((1i*w*T - 1)*exp(1i*w*T) + 1)/(1i*w*T^2)...
    - damp_coeff * P*exp(-damp_coeff*x)/(1i*w*T) * (exp(1i*w*T-1));
grad_v_in = P*((1i*w*T - 1)*exp(1i*w*T) + 1)/(1i*w*T^2*rho_l*c)...
    - damp_coeff * P*exp(-damp_coeff*x)/(1i*w*T*rho_l*c) * (exp(1i*w*T-1));

Fac = 4*pi/r * r^3 * (f_1 * 1/(2*rho_l*c^2) * grad_p_in - f_2 * 3/4*rho_l*grad_v_in);


end