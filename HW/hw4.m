%% Setup

%{ 
IMPORTANT: some variables are shared between phases, so
each section must be run selectively: 1->3 or 2->3, 3->1 or 4->1, etc.
this is a very disgusting implementation so... yeah
Improvement would be to have conditional statements that automatically
determine which scenario to use: this will loop through the script,
and at the end of each loop, it set {V_o_a, I_L_a, t_A} = {V_o_c,I_L_c, t_C}
respectively and iterate.
I shouldve tried to imitate what the RR_Ex10_20_boost_converter does but
that is too technical
%}
clear all, clc, close all
syms s t L C R V_i V_m V_o V_d I_L I_R I_C I_R
syms I_L_a I_L_b I_L_c V_o_a V_o_b V_o_c
disp("Hi i suck at this");
vars = [V_o,V_m,I_L,I_R,I_C];

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


% ON: t in [t_A,t_B].       V_o_a = V_o(t_C or t_A) , same for I_L_a
% OFF: t in [t_B, T_C].     V_o_b = V_o(t_B)        , same for I_L_b

f = 1.6e6; % Hz
D = 7/12;   % duty cycle with 1 being always ON
t_C = 1/f;
t_A = 0; t_B =(t_C-t_A)*D+t_A;
% generate times
tvalON = linspace(t_A,t_B,200); tON = [t_B,t_A];
tvalOFF= linspace(t_B,t_C,200); tOFF= [t_C,t_B];
% Init cond
subvars =   [V_i,V_d ,L    , C     ,R   ,V_o_a, I_L_a,V_o_b,I_L_b,V_o_c,I_L_c];
subvarsic = [5  ,.5  ,10e-6, 4.7e-6, 250, 0   , 0    ,0    ,0    , 0   , 0];
hold on;
%% phase 1 run: ON
A_1 = solve(eqns1,vars);
V_1 = ilaplace(A_1.V_o); % simple enough to do this
I_L1 = ilaplace(A_1.I_L);
% modify so that it has the time delay
V_1 = V_o_a*exp(-(t-t_A)/(C*R))
I_L1 = I_L_a + V_i/L*(t-t_A)
%insert variable values - FORMATTING: V_1 is with syms, V1 is with values
V1 = subs(V_1,subvars,subvarsic);
IL1 = subs(I_L1,subvars,subvarsic);
Vob1 = subs(V1,t,t_B) % generate IC's for next phase
ILb1 = subs(IL1,t,t_B)
subvarsic(8)= Vob1; % store in the sub vars IC array
subvarsic(9)= ILb1;
V1 = subs(V1,t,tvalON); % generate data points for duration of phase
hold on;
plot(tvalON,V1);
%% phase 2 run: ON
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
B0 = (b1-b2*sigma2)/omega2 - B2*sigma2/omega2

%I_L2: den (and hence a's) are the same, but num is different:
num2=numden(A_2.I_L); n2 = num2cell(simplify(flip(coeffs(num2,s))/(C*L*R)));
[c2 c1 c0] = deal(n2{:})
C2 = c0/(sigma2^2+omega2^2)
C1 = c2-C2
C0 = (c1-c2*sigma2)/omega2 - C2*sigma2/omega2
%{
V_0_2a = simplify(B2*1/s);
V_0_2b = simplify(B1*(s+sigma2)/((s+sigma2)^2+omega2^2));
V_0_2c = simplify(B0*omega2/((s+sigma2)^2+omega2^2));
V_0_2 = simplify(ilaplace(V_0_2a))+ simplify(ilaplace(V_0_2b))+ simplify(ilaplace(V_0_2c))
%}
I_2 = C2 + C1*exp(-sigma2*(t-t_A))*cos(omega2*(t-t_A))+C0*exp(-sigma2*(t-t_A))*sin(omega2*(t-t_A));
%V_o_2 = subs(V_o_2,subvars,subvarsic);
%v2 = subs(V_0_2,t,tvalON);
V_2 = B2 + B1*exp(-sigma2*(t-t_A))*cos(omega2*(t-t_A))+B0*exp(-sigma2*(t-t_A))*sin(omega2*(t-t_A));
V2 = subs(V_2,subvars,subvarsic); % substitude variable values
I2 = subs(I_2,subvars,subvarsic);
Vob2 = subs(V2,t,t_B) % solve for IC's at the end of the phase
ILb2 = subs(I2,t,t_B)

subvarsic(8)= Vob2; % store in the sub vars IC array
subvarsic(9)= ILb2;
V2 = subs(V2,t,tvalON); % generate data points for duration of phase

plot(tvalON,V2);

%% phase 3 run: OFF
A_3 = solve(eqns3,vars) %again very similar to 1
V_3 = V_o_b*exp(-(t-t_B)/(C*R));
% however, 
A_3.I_L
% returns 0 which is an issue for the setup of eqns3.
%% phase 4 run: OFF
A_4 = solve(eqns4,vars)
[num4,den4] = numden(A_4.V_o)
n4 = flip(coeffs(num4, s));   d4 = flip(coeffs(den4,s)); % extract coefficients in correct order
n4 = simplify(n4/d4(1));    d4 = simplify(d4/d4(1)); % do some maths
% warning: it drops the s*[rest of the stuff]
% in the form of V_o = [b2* s^2 + b1*s^1 + b0]/[s*(s^2 + a_1*s +a_o)]:
n4 = num2cell(n4); d4 = num2cell(d4);
[b2 b1 b0] = deal(n4{:});
[a2 a1 a0] = deal(d4{:});
sigma4 = a1/2;
omega4 = (a0-(a1^2)/4)^(1/2);
B2 = b0/(sigma4^2+omega4^2)
B1 = b2-B2
B0 = (b1-b2*sigma4)/omega4 - B2*sigma4/omega4

%I_L2: den (and hence a's) are the same, but num is different:
num4=numden(A_4.I_L); n4 = num2cell(simplify(flip(coeffs(num4,s))/(C*L*R)));
[c2 c1 c0] = deal(n4{:})
C2 = c0/(sigma4^2+omega4^2)
C1 = c2-C2
C0 = (c1-c2*sigma4)/omega4 - C2*sigma4/omega4

I_4 = C2 + C1*exp(-sigma4*(t-t_B))*cos(omega4*(t-t_B))+C0*exp(-sigma4*(t-t_B))*sin(omega4*(t-t_B));
%V_o_4 = subs(V_o_4,subvars,subvarsic);
%v4 = subs(V_0_4,t,tvalOFF);
V_4 = B2 + B1*exp(-sigma4*(t-t_B))*cos(omega4*(t-t_B))+B0*exp(-sigma4*(t-t_B))*sin(omega4*(t-t_B));
V4 = subs(V_4,subvars,subvarsic); % substitude variable values
I4 = subs(I_4,subvars,subvarsic);
Voc4 = subs(V4,t_C) % solve for IC's at the end of the phase
ILc4 = subs(I4,t_C)
subvarsic(6)= Voc4; % store in the sub vars IC array
subvarsic(7)= ILc4;
V4 = subs(V4,t,tvalOFF); % generate data points for duration of phase

plot(tvalOFF,V4);