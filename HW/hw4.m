syms s L C R V_i V_m V_o V_d I_L I_R I_C I_R

% start/phase A: t = [0, t_A], V_i = V_m, V_o = 0, button pressed.

eqn1a = V_m == V_i;         % V_o - V_m = 0 when there's no current through diode
eqn2 = V_m == L*s*I_L;      % V_m = L*d/dt(I_L)
eqn3 = I_C == C*s*V_o;      % I_C = C*d/dt(V_o-0)
eqn4 = V_o == 0;            % The resistor is isolated by the depletion zone in the diode and the capacitor
eqn5 = V_o == I_R*R;

% phase B: t = [t_A, t_B], button is released.
eqn1b = V_o - V_m == V_d;   % voltage drop from diode
eqn6 = I_L == I_C+I_R;      % KCL at V_o node
A_A = solve(eqn1a,eqn2,eqn2,eqn3,eqn4, eqn5, V_o, V_m, I_L, I_R, I_C);
A_A.V_o
A_B = solve(eqn1b,eqn2,eqn2,eqn3,eqn5,eqn6, V_o, V_m, I_L, I_R, I_C);
A_B.V_o
