% RR_Ex10_02_passive_filters
%% Q1
syms s I_L I_C I_R V_o V_i R_load L C; 
%x = {I_L, I_C, I_R, V_o}
A = [1  -1  -1  0;      % I_L - I_C - I_R = 0;
    L*s  0   0  1;      % L*s*I_L + V_o = V_i;
    0   -1   0  C*s;    % -I_C + C*s*V_o = 0;
    0    0   R_load  -1]; % I_R*R_load -V_o = 0;
b = [0; V_i; 0; 0];
x = A\b;
x=simplify(x);
V_o_undamped = x(4);
%% Q2
omega4 = 10;
fignum = 4;  % in the case i need to append this to RR_Ex10_02_passive_filters.m
legentries = {}; % to store each legend entry
for zeta = [0.1 0.7 1]
    F_undamped = RR_tf([omega4^2],[1 2*omega4*zeta omega4^2]);
    RR_bode(F_undamped);
    legentries{end+1} = ['\zeta = ' num2str(zeta)];
    legend(legentries);
    pause;
end

%% Q3
syms s I_L I_Cf I_Cd I_Rd V_o V_mid L Cd Cf Rd V_i
% x = {I_L, I_Cf, I_Cd, I_Rd, V_o, V_mid}
A= [1   -1  -1   0   0   0;     % I_L - I_Cf - I_Rd = 0
    0    0   1  -1   0   0;     % I_Cd - I_Rd = 0
    L*s  0   0   0   1   0;     % L*s*I_L + V_o = V_i
    0   -1   0   0  Cf*s 0;     % -I_Cf + Cf*s*V_o = 0
    0    0   0  Rd  -1   1;     % Rd*I_Rd - V_o + V_mid = 0
    0    0  -1   0   0  Cd*s];  % -I_Cd + Cd*s*V_mid = 0
b = [0;0;V_i;0;0;0];
x=A\b; V_o_PDF =simplify(x(5))

% Q4
Rd = sqrt(L/C); Cd = 4*Cf
F_PDF = RR_tf([Cd*Rd 1],[Cd*Cf*L*Rd (Cd+Cf)*L Cd*Rd 1]);
RR_bode(F_PDF);
%% In class example
%{
syms zeta
G1 = RR_tf([1 10 10000],[100 ,200* zeta + 1, 2*zeta+100, 1])
zeta = 0.1
G1 = RR_tf([1 10 10000],[100 ,200* zeta + 1, 2*zeta+100, 1])
RR_bode(G1);
%}