function [xe,A_d,B_d1,B_d2,C,i_b,Ge] = linear_matrix(p,Ts)

G_b = p(1);
V_I = p(2);
S_I = p(3);
k_tau = p(4);
k_d = p(7);
k_cl = p(8);
S_g = p(9);
V_g = p(10);
p_2 = p(11);
BW = p(12);
f_c = p(13);
I_b = p(14);

i_b = I_b*k_cl*V_I*BW;

Xe = S_I*(i_b/(k_cl*BW*V_I) - I_b);
Ge = (S_g*G_b)/(S_g+Xe);


A = ...
    [-(S_g+Xe), -Ge, 0, 0, 0;
    0, -p_2, 0, 0, (p_2*S_I)/(BW*V_I);
    0, 0, -k_d, 0, 0;
    0, 0, k_d, -k_d, 0;
    0, 0, 0, k_d, -k_cl
    ];
B_1 = [0; 0; 1; 0; 0];
B_2 = [1/V_g; 0; 0; 0; 0];
C = [1 0 0 0 0];


xe = [Ge; Xe; i_b/k_d; i_b/k_d; i_b/k_cl];


A_d = expm(A*Ts);
fun = @(tau) expm(A*tau);
B_d1 = integral(fun,0,Ts,"ArrayValued",true)*B_1;
B_d2 = integral(fun,0,Ts,"ArrayValued",true)*B_2;
end