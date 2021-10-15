function Fac = GorkovScript(c_l, rho_l, r, f_1, f_2, t, x)

j = sqrt(-1);

%The Mean Square of P and v
x0 = 5*10^-2; %5cm
T = 10; %10 sec
p_ms = 1/T * integral(@(t) (max(P(:,x0)) * exp(j*w*t-j*x0*2*pi/lambda - damp_coeff * x0)).^2,0,T);
v_ms = 1/T * integral(@(t) (1/(rho_l*c_l) * max(P(:,x0)) * exp(j*w*t-j*x0*2*pi/lambda - damp_coeff * x0)).^2,0,T);


Fac = 4*pi/3*r^3*(f_1*1/(2*rho_l*c_l^2)*p_ms-f_2*3/4*rho_l*v_ms);
end