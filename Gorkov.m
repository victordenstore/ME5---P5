
         function [v,p_mean_square,v_in,v_mean_square,Uac,Fac] = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B)
         

    syms x t 
    
P_field = P_surface*(exp(1i*w*t+1i*x*2*pi/lambda) + R1 * exp(1i*w*t-1i*x*2*pi/lambda))...
         + P_surface*(R2 * exp(1i*w*t+1i*x*2*pi/lambda) + exp(1i*w*t-1i*x*2*pi/lambda))*exp(1i*pi);

p_mean_square = 1/T * int(real(P_field)^2,t,[0 T]);

v_in = -1i/(w*rho_l)*diff(P_field,x);
v_mean_square = 1/T*int(real(v_in)^2,t,[0 T]);


Uac = 4*pi/3 * r^3 * (f_1 * 1/(2*rho_l*c^2) * p_mean_square - f_2 * 3/4*rho_l*v_mean_square);
Fac= -(diff(Uac,x));
v = Fac/B;
end