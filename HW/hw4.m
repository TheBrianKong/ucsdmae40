clc, clear all, close all
hold off;
syms s t L C R V_i V_m V_o V_d I_L I_R I_C I_R
syms I_L_a I_L_b V_o_a V_o_b
disp("Hi i suck at this");
vars = [V_o,V_m,I_L,I_R,I_C];
%% Sys of eqns
%1: phase A with no diode
eqn11 = V_m == V_i/s;
eqn12 = V_m == L*(s*I_L - I_L_a);
eqn13 = I_C == C*(s*V_o - V_o_a);
eqn14 = V_o == I_R*R;
eqn15 = I_C + I_R == 0;
eqns1 = [eqn11,eqn12,eqn13,eqn14,eqn15];
%2: phase A with diode
eqn21 = V_o - V_m == V_d/s;
eqn22 = V_m - 0 == L*(s*I_L - I_L_a);
eqn23 = I_C == C*(s*V_o - V_o_a);
eqn24 = -I_L == I_C+ I_R;
eqn25 = V_o == I_R*R;
eqns2 = [eqn21, eqn22,eqn23,eqn24,eqn25];
%3: phase B with no diode
eqn31 = V_m-0 == L*(s*I_L-I_L_b);
eqn32 = I_C+I_R == 0;       % V_o KCL
eqn33 = V_o == I_R*R;
eqn34 = -I_L == I_C + I_R;  % gnd KCL
eqn35 = I_C == C*(s*V_o - V_o_b);
eqns3 = [eqn31,eqn32,eqn33,eqn34,eqn35];
%4: phase B with diode
eqn41 = V_o - V_m == V_d/s;
eqn42 = V_m - 0 == L*(s*I_L - I_L_b);
eqn43 = I_C == C*(s*V_o - V_o_b);
eqn44 = -I_L == I_C + I_R;
eqn45 = V_o == I_R*R;
eqns4 = [eqn41,eqn42,eqn43,eqn44,eqn45];
%% Calc for phases

% ON: t in [t_A,t_B].       V_o_a = V_o(t_C or t_A) , same for I_L_a
% OFF: t in [t_B, T_C].     V_o_b = V_o(t_B)        , same for I_L_b
t_A = 0; t_B = 6e-3; t_C = 2e-3;
tvalON = linspace(t_A,t_B,200); tON = [t_B,t_A];
tvalOFF= linspace(t_B,t_C,200); tOFF= [t_C,t_B];
% Init cond
subvars =   [V_i,V_d ,L    , C     ,R   ,V_o_a, I_L_a];
subvarsic = [5  ,.5  ,10e-6, 4.7e-6, 250, 0   , 0];

% phase 1 run: ON
A_1 = solve(eqns1,vars);
V_1 = ilaplace(A_1.V_o); % simple enough to do this
I_L1 = ilaplace(A_1.I_L);
% modify so that it has the time delay
V_1 = V_o_a*exp(-(t-t_A)/(C*R))
I_L1 = I_L_a + V_i/L*(t-t_A)
%insert variable values - FORMATTING: V_1 is with syms, V1 is with values
V1 = subs(V_1,subvars,subvarsic);
IL1 = subs(I_L1,subvars,subvarsic);
V0b = subs(V_1,[t,t_A],tON) 
ILb = subs(I_L1,[t,t_A],tON) 
% phase 2 run: ON
A_2 = solve(eqns2,vars)
[num2,den2] = numden(A_2.V_o)
n2 = flip(coeffs(num2, s));   d2 = flip(coeffs(den2,s)); % extract coefficients in correct order
n2 = simplify(n2/d2(1));    d2 = simplify(d2/d2(1)); % do some maths
% warning: it drops the s*[rest of the stuff]
% in the form of V_o = [b2* s^2 + b1*s^1 + b0]/[s*(s^2 + a_1*s +a_o)]:
n2 = num2cell(n2); d2 = num2cell(d2);
[b2 b1 b0] = deal(n2{:});
[a2 a1 a0] = deal(d2{:});
sigma2 = a1/2;
omega2 = (a0-(a1^2)/4)^(1/2);
B2 = b0/(sigma2^2+omega2^2)
B1 = b2-B2
B0 = (b1-b2*sigma2)/omega2
%{
V_0_2a = simplify(B2*1/s);
V_0_2b = simplify(B1*(s+sigma2)/((s+sigma2)^2+omega2^2));
V_0_2c = simplify(B0*omega2/((s+sigma2)^2+omega2^2));
V_0_2 = simplify(ilaplace(V_0_2a))+ simplify(ilaplace(V_0_2b))+ simplify(ilaplace(V_0_2c))
%}

%V_0_2 = subs(V_0_2,subvars,subvarsic);
%v2 = subs(V_0_2,t,tval);
V_2 = B2 + B1*exp(-sigma2*(t-t_A))*cos(omega2*(t-t_A))+B0*exp(-sigma2*(t-t_A))*sin(omega2*(t-t_A));
V2 = subs(V_2,subvars,subvarsic);
V2 = subs(V2,t,tval);
hold on;
plot(tval,V2,'b');

% phase 3 run: OFF
A_3 = solve(eqns3,vars)
% phase 4 run: OFF
A_4 = solve(eqns4,vars)
[num4,den4] = numden(A_4.V_o)
n4 = flip(coeffs(num4, s));   d4 = flip(coeffs(den4,s)); % extract coefficients in correct order
n4 = simplify(n4/d4(1));    d4 = simplify(d4/d4(1)); % do some maths
% warning: it drops the s*[rest of the stuff]
% in the form of V_o = [b2* s^2 + b1*s^1 + b0]/[s*(s^2 + a_1*s +a_o)]:
n4 = num2cell(n4); d4 = num2cell(d4);
[b2 b1 b0] = deal(n2{:});
[a2 a1 a0] = deal(d2{:});
sigma4 = a1/2;
omega4 = (a0-(a1^2)/4)^(1/2);
B2 = b0/(sigma4^2+omega4^2)
B1 = b2-B2
B0 = (b1-b2*sigma4)/omega4
%{
V_0_4a = B2*1/s;
V_0_4b = B1*(s+omega4)/((s+sigma4)^2+omega4^2)
V_0_4c = B0*omega4/((s+sigma4)^2+omega4^2)
simplify(ilaplace(V_0_4a))
simplify(ilaplace(V_0_4b))
simplify(ilaplace(V_0_4c))
%}