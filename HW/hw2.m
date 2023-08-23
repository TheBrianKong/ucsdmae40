%% Q1
syms s V_1 V_2 I_o I_a I_b I_c I_d R_a C_b C_c R_d V_o
%x={V_1 V_2 Io    Ia         Ib     Ic      Id}
A = [1  0   0      R_a       0      0       0 ; % V1 + Ia*Ra = Vo
     0  1   0      0      1/(C_b*s) 0       0 ; % V2 - Vo = Ib/(C_b*s ) <- % Cb*s * V2 + Ib = Cb*s* Vo
     1 -1   0      0         0   -1/(C_c*s) 0 ; % Cc*s * V1 - Cc*s* V2 - Ic = 0
     0  1   0      0         0      0       -R_d; % V2 - Rd * Id = 0
     0  0   1     -1        -1      0       0 ; % Io - Ia - Ib = 0
     0  0   0      1         0     -1       0 ; % Ia - Ic = 0
     0  0   0      0         1      1      -1 % Ib + Ic - Id = 0
    ];
b = [V_o; V_o; 0;0;0;0;0];
x= A\b;
F1 = simplify(x(1))
% latex(x(1)) % for personal use to check when formatted.
syms R C; % Q1c: R_a = 2R, R_d = R/2, C_b = C_c = C
F1=subs(F1,{R_a, R_d, C_b, C_c},{2*R, R/2, C, C})
% F1 = V_1/V_o = funct(C,R,s)
poles(F1);
%F1 = RR_tf(F1)
% From prior calculations: V_1/V_o = (C^2*R^2*s^2 + C*R*s + 1))/(C^2*R^2*s^2 + 3*C*R*s + 1)
F1_num = [C^2*R^2 C*R 1]; F1_den = [C^2*R^2 3*C*R 1];
% which turns into the following:
omega1 = 10;
hold on; grid on;
F1_num = [1 omega1 omega1^2]; F1_den = [1 3*omega1 omega1^2];
bode(tf(F1_num,F1_den)); 
legentries = {['\zeta = 3/2']}; % to store each legend entry
for zeta = 10.^[-2:1]
F1_num = [1 omega1 omega1^2]; F1_den = [1 2*zeta*omega1 omega1^2];
bode(tf(F1_num,F1_den));
    legentries{end+1} = ['\zeta = ' num2str(zeta)];
    legend(legentries,'Location','se');
end
%% Q2
syms s V_1 V_2 I_0 I_a I_b I_c I_d C_a R_b R_c C_d V_0
%x= [V_1 V_2    I_0 I_a        I_b  I_c     I_d]
A = [0   0      1   -1         -1   0       0       % I_0-Ia-Ib=0
     0   0      0   1          0    -1      0       % Ia -Ic = 0
     0   0      1   0          0    0       -1      % Ib + Ic = Id = I_0
     1   0      0   1/(C_a*s)  0    0       0       % V_0-V_1 = Ia/(Ca*s)
     0   1      0   0          R_b  0       0       % V_2 + Ib*Rb = V_0
     1  -1      0   0          0    -R_c    0       % V_1 - V_2 -Ic*Rc = 0
     0  C_d*s   0   0          0    0       -1      % Cd*s*V2 - Id = 0
    ];
b=[0;0;0;V_0;V_0;0;0];
x = A\b; % expand((R_c*C_a*s*(V_1-V_0)+V_1)*(C_d*s+1/R_b+1/R_c)) % grunt work
F2 = simplify(x(1))
syms C R
F2 = simplify(subs(F2,{C_a,R_b,R_c,C_d},{C/2,R,R,2*C}))
F1_num = [C^2*R^2 C*R 1]; F1_den = [C^2*R^2 3*C*R 1];
% which turns into the following:
omega2 = 10;
hold on; grid on;
F2_num = [1 omega2 omega2^2]; F2_den = [1 3*omega2 omega2^2];
bode(tf(F2_num,F2_den));
legentries = {['\zeta = 3/2']}; % to store each legend entry
for zeta = 10.^[-2:1]
bode(tf([1 omega2 omega2^2],[1 2*zeta*omega2 omega2^2]));
    legentries{end+1} = ['\zeta = ' num2str(zeta)];
    legend(legentries,'Location','se');
end