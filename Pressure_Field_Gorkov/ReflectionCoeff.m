syms n
P=symsum((0.9158^(2*n)*exp(-2*n*0.00448*0.3833)),n,1,inf)
syms n
P_1=symsum((0.9158^(2*n-1)*exp(-2*(n-1)*0.00448*0.3833)),n,1,inf)
