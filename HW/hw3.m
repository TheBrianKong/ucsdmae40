%% Q1
clear, clc, close all
syms s L V_m C V_o R I_L I_C I_R I_load c_1 V_i

%x= {V_o    V_m     I_L    I_C     I_R     I_load}

A = [0      1       L*s     0       0       0;      % L*s*I_L + V_m = Vi
     -1     1       0    -1/(C*s)   0       0;      % I_C = C*s(V_m - V_o)
     1      0       0       0       -R      0;      % V_o = I_R*R
     -1     0       0       0       0       R/c_1   % V_o = I_load*R/c_1
     0      0       1       -1      0       0;      % I_L - I_C = 0
     0      0       0       -1      1       1;      % I_C = I_R + I_load
    ];
b = [V_i; 0; 0; 0; 0; 0];

x = A\b;
V_o_resp = simplify(x(1))
%% Optional solve for terms
%{
Here we get that V_o/V_i = [C*R*V*s] / [C*L*(c_1+1)*s^2 + C*R*s + (1+c_1)]
We can simplify this with the following:
%}
[num, den]= numden(V_o_resp/V_i); % divide by V_i and extract num and den
num = flip(coeffs(num, s)); den = flip(coeffs(den,s)); % extract coefficients in correct order
num = simplify(num/den(1)); den = simplify(den/den(1));% clean up expressions
omega_0 = sqrt(den(3))         % extract values of omega_0 and Q from the form:
Q = simplify(omega_0/den(2))    % den =  [1 omega_0/Q omega_0^2]
% further simplification yields: Q = (c_1+1)/R*(L/C)^(1/2)
%% Plotting
syms omega_0 Q
num1= [omega_0/Q 0]; den1 = [1 omega_0/Q omega_0^2];
% F_1 = RR_tf(num1, den1);
% substitute values for Q=5 and omega_0=10
F_1 = RR_tf(subs(num1,[Q omega_0],[5 10]), ...
            subs(den1,[Q omega_0],[5 10]));
RR_bode(F_1)
