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
dT_pure = matlabFunction(dTi);

Tarr = [T1,T2,T3,T4,T5,T6]
Tnames = ['T_1','T_2','T_3','T_4','T_5','T_6']

nT = max(size(Tarr))
nCells = nT-1 % We have one more T (surface) than cells


u_params = [U_A] % Note! You will have do do manual changes if you increase params here by adding more u_UA etc functions
fluctating_params = [Qdotger, mdot,Delta_P_H2O,T_ar]

A = jacobian(dTi, Tarr)
%B_full = jacobian(dTi, [Qdotger,U_A,mdot,Delta_P_H2O, T_ar])
B = jacobian(dTi, u_params)
E = jacobian(dTi, fluctating_params)

latexify("\frac{d}{dt}\vec T =",dTi)
latexify("A=",A)
latexify("B=",B)

% y independent of x after a substitution, i.e. we cant update the value again

LiterToM3 = 1000
F = 0.9
cH2O = 4186 % J/kgK
rho = 1*LiterToM3 % kg/L
mdot = 25/60 % kg/s
Qdotger = 6000 % W

Delta_P_H2O = 35000 % Pa
Delta_P_R = Delta_P_H2O/nCells % Pa

V_motor = 1.5/LiterToM3 % L
V_cell = V_motor/nCells

F = 0.9
T_ar = 298 % K
U_A = 100/nCells % W/K

Qdotger_fluc = 1500 % W
mdot_fluc = 5/60 % kg/s
Delta_P_H20_fluc = 5000 % Pa
T_ar_fluc = 10 % C

fluc_arr = [Qdotger_fluc, mdot_fluc, Delta_P_H20_fluc, T_ar_fluc]

vd = diag(fluc_arr./fluctating_params) % really important that these are in right order

vd = subs(vd)
dTi = subs(dTi);

% Equilibrium conditions

dT_handle = matlabFunction(dTi);
dT_lambda = @(t,T)dT_handle(T(1),T(2),T(3),T(4),T(5),T(6))
%[t,y] = ode45(dT_lambda(),[0 200],[400;399;398;397;396;395]);

tSimulation = 0:1:1000;
TarrSimulation = [370;369;368;367;366;365]


% TODO IMPLEMENT THE CONTROLLER LIKE THIS
%x0 = [-1; 0; pi+.1; 0];  % initial condition 
%wr = [1; 0; pi; 0];      % reference position
%u=@(x) U_A-K*(x - wr);       % control law
%u=@(x)-K*(x - wr);       % control law
%[t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);



eq_sol = solve(dTi)
T1 = eq_sol.T1;
T2 = eq_sol.T2;
T3 = eq_sol.T3;
T4 = eq_sol.T4;
T5 = eq_sol.T5;
T6 = eq_sol.T6;
T_equilibrium = [T1;T2;T3;T4;T5;T6]

A = double(subs(A));
B = double(subs(B));
E = double(subs(E));


C_full_observe = eye(size(A)); % Matriz de identidade
C = zeros([2, nT]); % We only measure 2 values
C(1) = 1; % We can only measure edge temperatures
C(end) = C(1);
%C(2,2 ) = 1

D_full_observe = zeros([size(C_full_observe,1),size(B,2)]); % Matriz de zeros
D = zeros([size(C,1),size(B,2)]); % Matriz de zeros


% Definindo o sistema
sys = ss(A, B, C_full_observe, D_full_observe)

%%

% Plot a step resonse 
t = 0:0.1:20;

respOpt = RespConfig;
respOpt.InputOffset = [0];
respOpt.Amplitude = [-2];
respOpt.InitialState = [0,0,0,0,0,0];
respOpt.Delay = 0;

figure;
step(sys,t,respOpt)
title('Stepresponse to changed U_A, open loop');

%%
% Obtendo os polos
polos = eig(A);
disp('Polos:');
disp(polos);

% Loop para calcular as funcoes de transferencia
for IU = 1:size(B,2)
    [num,den] = ss2tf(A, B, C, D, IU);
    num = num(end,:); % Converte a matriz num em um vetor de linha
    den = den(end,:); % Converte a matriz den em um vetor de linha
    TF = tf(num,den);
    disp(['Função de Transferência para a entrada ', num2str(IU), ':']);
    disp(TF);
end

% Plotando o diagrama de Bode
%figure;
%bode(sys);
%title('Diagrama de Bode');

%%

%%%% REGULADOR LINEAR QUADRÁTICO %%%%

% Definindo as matrizes de peso Q e R
Q = eye(size(A)); % Matriz de peso para os estados
R = 0.0001; % Peso para a entrada, %% Make it expensive to change fan speed

% Calculando o ganho do controlador LQR
Klq = lqr(A, B, Q, R);

% Sistema em malha fechada com controlador LQR
A_cl1 = A - B * Klq;
sys_cl1 = ss(A_cl1, B, C_full_observe, D_full_observe);

% Exibindo os polos
polos_lqr = eig(A_cl1);
disp('Polos do sistema com controlador LQR:');
disp(polos_lqr);

t = 0:0.1:20;
figure;
step(sys_cl1,t,respOpt)
title('Stepresponse to changed U_A, LQR');
ylabel('\Delta T compared to equlibrium (K)');

%% Pole Placement

% Uses LQR poles and tries to vary them slightly to manually
% decide if there are more suitable placements

%polesFactorArr = [0.3, 0.8, 0.9,1, 1.1, 1.2, 3]
polesFactorArr = [0.3]


for poleFactor = 1:size(polesFactorArr,2)
    desired_poles = polos_lqr * polesFactorArr(poleFactor)
    Kp = place(A, B, desired_poles);
    A_poleTemp = A - B * Kp;
    disp(eig(A_poleTemp))
    sys_poleTemp = ss(A_poleTemp, B, C_full_observe, D_full_observe);
    % Add plotting code here @ Arthur, plot the U
    figure;
    step(sys_poleTemp,t,respOpt)
    title(append('Step response to changed U_A, Pole placement: LQR \times' ,num2str(polesFactorArr(poleFactor))));
    ylabel('\Delta T compared to equlibrium (K)');
end

%% LQE



sys_cl_k = ss(A, B, C, D);
sys_cl_k.InputName = 'U_A';
sys_cl_k.OutputName = 'T1-6';

%[L,P,E] = lqe(A,Vd,C,Vd,Vn)
%Kf = (lqr(A',C',Vd,Vn))'

Q = 10;
R = 10;
N = 0; % No correlation

[kalmf,L,P] = kalman(sys_cl_k,Q,R,N);
size(kalmf)

kalmf.InputGroup

figure;
disp("not using kalman filter")
step(sys_cl_k,t,respOpt)
title('Stepresponse to changed U_A, Kalman');


% TODO add simulation with kalman filter, basically what we have here
% https://www.mathworks.com/help/control/ug/kalman-filtering.html
% Also, test controlability
% Pole placement just change LQR with 10%
% Plot all poles




%% Simulations performed here tambem


u_UA_malha_abert = @(T)U_A%-K*(T-Tarr)
[t,y_malha_abert] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_malha_abert(T),V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);

u_UA_LQR = @(T)double(-Klq*(T-T_equilibrium))
[t,y_LQR] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_LQR(T),V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);


figure;
plot(tSimulation,y_malha_abert,'r')
hold on
plot(tSimulation,y_LQR,'b')
hold on











%% Simulation part IGNORE THIS


%{
% Simular a resposta ao degrau do sistema
[ylq, t, xlq] = lsim(sys_cl1, u, t, x0);

% Plotar a resposta do sistema
figure;
plot(t, ylq)
xlabel('Tempo (s)')
ylabel('Saída do sistema (K)')
legend(Tnames)
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

%%%% MALHA FECHADA %%%%

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
%}
%%%% ALOCACAO DE POLOS %%%%
%{
% Polos desejados
desired_poles = [-1 + 1i, -1 - 1i , -0.0001];
%TODO update based on LQR-search?


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



% Plotando os polos
figure;
pzmap(sys,'r',sys_cl, 'g',sys_cl1,'b');
legend("Sem alocação", "Com alocação", "LQR")

title('Mapa de Polos e Zeros');
%}