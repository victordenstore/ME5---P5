         function Fac = Gorkov(c, rho_l, r, f_1, f_2,w,t,x,damp_coeff,P,lambda)

%The Mean Square of P and v

% T = 2*pi/w; % one period
% p_ms = 1/T * integral(@(t) (P * exp(1i*w*t-1i*x*2*pi/lambda - damp_coeff * x).^2),0,T);
% v_ms = 1/T * integral(@(t) (1/(rho_l*c) * P * exp(1i*w*t-1i*x*2*pi/lambda - damp_coeff * x).^2),0,T);
%  
% 
% Uac = 4*pi/3*r^3*(f_1*1/(2*rho_l*c^2)*p_ms-f_2*3/4*rho_l*v_ms);
% j = sqrt(-1);
% num1 = 0.25 * j * f_1 * P^2 * (exp(2*j*w*t - 8.377580408 * j *x) - exp(-8.377580408*j*x)) * exp(-2 * damp_coeff * x);
% denom1 = rho_l * c^2 * t^2 * w;
% num2 = 0.5 * f_1 * P^2 * exp(2*j*w*t-8.377580408*j*x) * exp(-2*damp_coeff*x);
% denom2 = rho_l * c^2 * t;
% num3 = 0.375 * j * f_2 * P^2 *(-exp(2*j*w*t-8.377580408*j*x)+exp(-8.377580408*j*x))*exp(-2*damp_coeff*x);
% denom3 = rho_l * t^2 * c^2 * w;
% num4 = 0.75 * f_2 * P^2*exp(2*j*w*t-8.377580408*j*x)*exp(-2*damp_coeff*x);
% denom4 = rho_l * t * c*2;
% Fac = 1/3 * (4*pi*r^3*(num1/denom1 + num2/denom2 + num3/denom3 - num4/denom4));

p_ms = P * exp(1i*w*t-1i*x*2*pi/lambda - damp_coeff * x);
v_ms = 1/(rho_l*c) * P * exp(1i*w*t-1i*x*2*pi/lambda - damp_coeff * x);
Fac = 4*pi/3*r^3*(f_1*1/(2*rho_l*c^2)*abs(p_ms^2)-f_2*3/4*rho_l*abs(v_ms)^2);  %% Not sure about p_ms and v_ms. And imag part of f1 and f2?


end