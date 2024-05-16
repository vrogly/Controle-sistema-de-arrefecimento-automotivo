clc;clear;
close all

disp("Running")

syms T1 T2 T3 T4 T5 T6 T_ar mdot cH2O Qdotger rho Delta_P_R Delta_P_H2O V_motor V_cell U_A F

dTi = [(mdot*cH2O*(T6-T1) + Qdotger +(mdot /rho)*Delta_P_H2O ) / (rho * V_motor * cH2O), %
    (-F*U_A*((T1-T2)/(log((T1-T_ar)/(T2-T_ar)))) + mdot*cH2O*(T1-T2) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),%
    (-F*U_A*((T2-T3)/(log((T2-T_ar)/(T3-T_ar)))) + mdot*cH2O*(T2-T3) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T3-T4)/(log((T3-T_ar)/(T4-T_ar)))) + mdot*cH2O*(T3-T4) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T4-T5)/(log((T4-T_ar)/(T5-T_ar)))) + mdot*cH2O*(T4-T5) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T5-T6)/(log((T5-T_ar)/(T6-T_ar)))) + mdot*cH2O*(T5-T6) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    ]


Tarr = [T1,T2,T3,T4,T5,T6]
nCells = max(size(Tcells))-1 % We have one more T (surface) than cells


A = jacobian(dTi, Tarr)
B_full = jacobian(dTi, [Qdotger,U_A,mdot,Delta_P_H2O, T_ar])
B = jacobian(dTi, [mdot, U_A])

latexify("\frac{d}{dt}\vec T =",dTi)
latexify("A=",A)
latexify("B=",B)

% y independent of x after a substitution, i.e. we cant update the value again

LiterToM3 = 1000
F = 0.9
cH2O = 4186 % J/kgK
rho = 1*LiterToM3 % kg/L
mdot = 25/60 % kg/s
Qdotger = 7500 % W

Delta_P_H2O = 35000 % Pa
Delta_P_R = Delta_P_H2O/nCells % Pa

V_motor = 1.5/LiterToM3 % L
V_cell = V_motor/nCells

%V_I = 0.75/LiterToM3 % L
%V_II = 0.75/LiterToM3 % L
F = 0.9
T_ar = 298 % K
U_A = 80 % W/K

dTi = subs(dTi)

% Equilibrium conditions

eq_sol = solve(dTi)
T1 = eq_sol.T1;
T2 = eq_sol.T2;
T3 = eq_sol.T3;
T4 = eq_sol.T4;
T5 = eq_sol.T5;
T6 = eq_sol.T6;


A = subs(A)
B_full = subs(B_full)
B = subs(B)

C_full_observe = eye(size(A)); % Matriz de identidade
C = [[1,0,0],[0,0,0],[0,0,1]]

D1 = zeros(size(C_full_observe,1), size(B_full,2)); % Matriz de zeros

% Definindo o sistema
sys = ss(A, B, C1_full_observe, D1);

% Obtendo os polos
polos = eig(A);
disp('Polos:');
disp(polos);

% Loop para calcular as funcoes de transferencia
for IU = 1:size(B1,2)
    [num,den] = ss2tf(A, B1, C1, D1, IU);
    num = num(end,:); % Converte a matriz num em um vetor de linha
    den = den(end,:); % Converte a matriz den em um vetor de linha
    TF = tf(num,den);
    disp(['Função de Transferência para a entrada ', num2str(IU), ':']);
    disp(TF);
end

% Plotando o diagrama de Bode
figure;
bode(sys);
title('Diagrama de Bode');



%%%% MALHA FECHADA %%%%

B = [0; e_dUA_term; f_dUA_term];
disp(B); % Matriz de controle

C = eye(size(A)); % Matriz de identidade
D = zeros(size(A, 1), size(B, 2)); % Matriz de zeros

E = [d_dQger, d_dm_dot, d_dDelta_P_H2O, 0;
     0, e_dm_dot, e_dDelta_P_H2O, e_dT_ar;
     0, f_dm_dot, f_dDelta_P_H2O, f_dT_ar];
disp(E); % Matriz de distúrbios


% Verificando a controlabilidade
Co = ctrb(sys);
rank_co = rank(Co);

if rank_co == size(A, 1)
    disp('O sistema é controlável.')
else
    disp('O sistema não é controlável.')
end

% Verificando a observabilidade
Ob = obsv(sys);
rank_ob = rank(Ob);

if rank_ob == size(A, 1)
    disp('O sistema é observável.')
else
    disp('O sistema não é observável.')
end

%%%% ALOCACAO DE POLOS %%%%

% Polos desejados
desired_poles = [-1 + 1i, -1 - 1i , -0.0001];

% Calcular o ganho do controlador K
Kp = place(A, B, desired_poles);
disp('A matriz de ganhos K é:');
disp(Kp);

% Novo sistema com realimentação de estado
A_cl = A - B * Kp;
sys_cl = ss(A_cl, B, C, D);

% Simular a resposta ao degrau do sistema
t = 0:0.01:15;
u = zeros(size(t));
x0 = [T1; T2; T3];
[yp, t, xp] = lsim(sys_cl, u, t, x0);

% Plotar a resposta do sistema
figure;
plot(t, yp)
xlabel('Tempo (s)')
ylabel('Saída do sistema (K)')
title('Resposta do Sistema com Alocação de Polos')
legend('T_1','T_2','T_3')
grid on

% Número de entradas do sistema
num_inputs = size(B1, 2);

% Criar um vetor de entrada com o mesmo número de colunas que o número de entradas do sistema
u1 = zeros(length(t), num_inputs);

% Agora você pode usar este vetor de entrada com a função lsim
[ynp, t, xnp] = lsim(sys, u1, t, x0);

% Plotar a resposta do sistema
figure;
plot(t, ynp)
xlabel('Tempo (s)')
ylabel('Saída do sistema (K)')
title('Resposta do Sistema sem Alocação de Polos')
legend('T_1','T_2','T_3')
grid on

%%%% REGULADOR LINEAR QUADRÁTICO %%%%

% Definindo as matrizes de peso Q e R
Q = eye(3); % Matriz de peso para os estados
R = 1; % Peso para a entrada

% Calculando o ganho do controlador LQR
Klq = lqr(A, B, Q, R);

% Sistema em malha fechada com controlador LQR
A_cl1 = A - B * Klq;
sys_cl1 = ss(A_cl1, B, C, D);

% Exibindo os polos
polos_lqr = eig(A_cl1);
disp('Polos do sistema com controlador LQR:');
disp(polos_lqr);

% Simular a resposta ao degrau do sistema
[ylq, t, xlq] = lsim(sys_cl1, u, t, x0);

% Plotar a resposta do sistema
figure;
plot(t, ylq)
xlabel('Tempo (s)')
ylabel('Saída do sistema (K)')
legend('T_1','T_2','T_3')
title('Resposta do Sistema com Controlador LQR')
grid on

% Comparacao entre LQR e Alocacao de polos 
figure;
plot(t, yp, '--', t, ylq);
xlabel('Tempo (s)');
ylabel('Saída do sistema (K)');
legend('Alocação de Polos T_1', 'Alocação de Polos T_2', 'Alocação de Polos T_3','LQR T_1','LQR T_2','LQR T_3');
title('Comparação entre Alocação de Polos e LQR');
grid on;

% Plotando os polos
figure;
pzmap(sys,'r',sys_cl, 'g',sys_cl1,'b');
legend("Sem alocação", "Com alocação", "LQR")

title('Mapa de Polos e Zeros');
