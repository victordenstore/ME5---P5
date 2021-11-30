% Within the mainfile of the model 'ModelScript', there is called upon this function to calculate various parameters such as velocity, force on partcile, etc.

         
function [v,p_mean_square,v_in,v_mean_square,Uac,Fac] = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B)
         

    syms x t                                                                        % In order to calculate the time-averaged pressure and the gradient of 
                                                                                    % Gor'kov's potential, x and t must be symbols which will later be replaced 
                                                                                    % by numbers for time and distance.
    
P_field = P_surface*(R1+1) * exp(1i*w*t-1i*x*2*pi/lambda) ...     % The standing pressure field wave expression
         - R2*P_surface * exp(1i*w*t+1i*(x)*2*pi/lambda);

p_mean_square = 1/T * int(real(P_field)^2,t,[0 T]);                                     % Time averaged pressure

v_in = -1i/(w*rho_l)*diff(P_field,x);
v_mean_square = 1/T*int(real(v_in)^2,t,[0 T]);


Uac = 4*pi/3 * r^3 * (f_1 * 1/(2*rho_l*c^2) * p_mean_square - f_2 * 3/4*rho_l*v_mean_square);   % Gor'kov's potential
Fac= -(diff(Uac,x));                                                                            % Force on air bubble
v = Fac/B;                                                                                      % velocity of air bubble
end