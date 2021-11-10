%% local comparison for meeting

x = input('put in whatever value you want here to compare: ');
[v,p_mean_square,v_in,v_mean_square,Uac,Fac]  = Gorkov(P_surface,w,lambda,R1,R2,T,rho_l,r,f_1,f_2,c,B);
Fac_local = double(subs(Fac));
v_local = double(subs(v));
clc
string1 = ['Fac_local = ' num2str(Fac_local)]
string2 = ['v_local = ' num2str(v_local)];
disp(string1); disp(string2);