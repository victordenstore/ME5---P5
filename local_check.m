%% local comparison for meeting

x = input('put in whatever value you want here to compare: ');
[v,p_mean_square,v_in,v_mean_square,Uac,Fac]  = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B);
Fac_local = double(subs(Fac));
v_local = double(subs(v));
clc
disp(Fac_local); disp(v_local);