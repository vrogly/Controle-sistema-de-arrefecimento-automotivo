clc;clear;
close all
set(groot,'defaultAxesFontSize',16)
set(groot,'DefaultFigureRenderer','painters');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

syms T1 T2 T3 T4 T5 T6 T_ar mdot cH2O Qdotger rho Delta_P_R Delta_P_H2O V_motor V_cell U_A F;

dTi = [(mdot*cH2O*(T6-T1) + Qdotger +(mdot /rho)*Delta_P_H2O ) / (rho * V_motor * cH2O), %
    (-F*U_A*((T1-T2)/(log((T1-T_ar)/(T2-T_ar)))) + mdot*cH2O*(T1-T2) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),%
    (-F*U_A*((T2-T3)/(log((T2-T_ar)/(T3-T_ar)))) + mdot*cH2O*(T2-T3) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T3-T4)/(log((T3-T_ar)/(T4-T_ar)))) + mdot*cH2O*(T3-T4) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T4-T5)/(log((T4-T_ar)/(T5-T_ar)))) + mdot*cH2O*(T4-T5) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    (-F*U_A*((T5-T6)/(log((T5-T_ar)/(T6-T_ar)))) + mdot*cH2O*(T5-T6) -(mdot/rho)*Delta_P_R ) / (rho * V_cell * cH2O),
    ];
dT_pure = matlabFunction(dTi);

Tarr = [T1,T2,T3,T4,T5,T6];
Tnames = ['T_1','T_2','T_3','T_4','T_5','T_6'];

nT = max(size(Tarr));
nCells = nT-1; % Temos um T (superfície) a mais que células


u_params = [U_A]; % Observação! Você terá que fazer alterações manuais se aumentar os parâmetros aqui adicionando mais funções u_UA etc.
fluctating_params = [Qdotger, mdot,Delta_P_H2O,T_ar];

A = jacobian(dTi, Tarr);
B = jacobian(dTi, u_params);
E = jacobian(dTi, fluctating_params);

latexify("\frac{d}{dt}\vec T =",dTi)
latexify("A=",A)
latexify("B=",B)

LiterToM3 = 1000;
F = 0.9;
cH2O = 4186; % J/kgK
rho = 1*LiterToM3; % kg/L
mdot = 25/60; % kg/s
Qdotger = 6000; % W

Delta_P_H2O = 35000; % Pa
Delta_P_R = Delta_P_H2O/nCells; % Pa

V_motor = 1.5/LiterToM3; % L
V_cell = V_motor/nCells;

F = 0.9; % -
T_ar = 298; % K
U_A = 100/nCells; % W/K 

dTi = subs(dTi);

% Simulação usada para comparação
dt = 0.01;
tSimulation = 0:dt:30;
TarrSimulation = [370;369;368;367;366;365];

eq_sol = solve(dTi);
T1 = eq_sol.T1;
T2 = eq_sol.T2;
T3 = eq_sol.T3;
T4 = eq_sol.T4;
T5 = eq_sol.T5;
T6 = eq_sol.T6;
T_equilibrium = [T1;T2;T3;T4;T5;T6];

A = double(subs(A));
B = double(subs(B));
E = double(subs(E));

C_full_observe = eye(size(A)); % Matriz de identidade

% Construa uma matriz onde observamos apenas as temperaturas selecionadas 
% O padrão é observar apenas a temperatura nas bordas
nMeasuresC = 1; 
C = zeros([nMeasuresC, nT]);
C(1) = 1;
%C(end) = C(1);
for nn = 2:nMeasuresC-1
    C(nn,nn) = 1;
end

D_full_observe = zeros([size(C_full_observe,1),size(B,2)]); % Matriz de zeros
D = zeros([size(C,1),size(B,2)]); % Matriz de zeros

s = tf('s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Controle moderno                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definindo o sistema malha aberta, observando todas as variáveis
sys_ol = ss(A, B, C_full_observe, D_full_observe)

% Obtendo os polos
polos = eig(A);
disp('Polos:');
disp(polos);

% Plotando o diagrama de Bode
figure;
bode(sys_ol);
title('Diagrama de Bode');


%%%% REGULADOR LINEAR QUADRÁTICO (LQR) %%%%

% Definindo as matrizes de peso Q e R
% Não adiantou colocar 0,1 ou 0,5 em vez de 0, só piorou
Q = diag([1,0,0,0,0,1]);%
% Torne caro alterar a velocidade do ventilador
R = 0.01; 

% Calculando o ganho do controlador LQR
Klqr = lqr(A, B, Q, R);
% Sistema em malha fechada com controlador LQR
A_cl_lqr = A - B * Klqr;

sys_cl_lqr = ss(A_cl_lqr, B, C_full_observe, D_full_observe);
Co = ctrb(sys_cl_lqr);
rank_co = rank(Co);

% Exibindo os polos
polos_lqr = eig(A_cl_lqr);
disp('Polos do sistema com controlador LQR:');
disp(polos_lqr);

% Compare malha aberta com malha fechada com LQR
u_UA_malha_aberta = @(T)U_A;
[t,y_malha_aberta] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_malha_aberta(T),V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);

u_UA_LQR = @(T)double(U_A-Klqr*(T-T_equilibrium));
[t,y_LQR] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_LQR(T),V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);

figure;
title("Controlando as condições iniciais para o estado desejado")
plot(tSimulation,y_malha_aberta,'r');
hold on;
plot(tSimulation,y_LQR,'b');
hold on;
legend("Malha aberta","","","","","","","LQR")
ylabel('$T_1, T_2, T_3, T_4, T_5, T_6$ (K)');
xlabel('Tempo (s)');
saveas(gcf,'aberta_lqr_T.eps','epsc')

u_LQR_history = [];
for i = 1:size(y_LQR,1)
    u_LQR_history = [u_LQR_history nCells * u_UA_LQR(transpose(y_LQR(i,:)))];
end    
figure;
title("Controlando as condições iniciais para o estado desejado")
plot(tSimulation, nCells*U_A*ones(size(tSimulation)),'r');hold on
plot(tSimulation,u_LQR_history,'b')
ylabel('$U_A$ (W/K)');
xlabel('Tempo (s)');
legend("Malha aberta","LQR");
saveas(gcf,'aberta_lqr_u.eps','epsc')


%%%% Alocação de Polos %%%%
%%
% Usa pólos LQR e tenta variá-los ligeiramente para decidir manualmente se há posicionamentos mais adequados
% Trace o controle ao longo do tempo para ver se ele não excede nenhuma especificação

polesFactorArr = [0.7, 0.9, 1, 1.1, 1.3];
semimanualPoleArr = [polos_lqr diag([1 1.01 1 1 1 1.01])*real(polos_lqr) diag([0.9 1 1 1 1 1])*polos_lqr diag([1.1 1 1 1 1 1])*polos_lqr diag([1.3 1 1 1 1 1])*polos_lqr]
T_pole_history = figure; 
u_pole_history = figure;

for poleFactor = 1:size(semimanualPoleArr,2)
    desired_poles = semimanualPoleArr(:,poleFactor)
    Kp = place(A, B, desired_poles);
    A_poleTemp = A - B * Kp;
    disp(eig(A_poleTemp))
    
    % Aqui poderíamos construir e simular um modelo de espaço de estados como
    % sys_poleTemp = ss(A_poleTemp, B, C_full_observe, D_full_observe);
    % Mas, em vez disso, simulamos dinâmicas não lineares completas e atualizamos o
    % variável de controle usando seu modelo controle linear

    u_UA_KP = @(T)double(U_A-Kp*(T-T_equilibrium));
    [t,y_KP] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger, ...
        T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_KP(T), ...
        V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);    
    
    figure(T_pole_history);
    plot(tSimulation,y_KP(:,1)); hold on;
    
    u_KP_history = [];
    for i = 1:size(y_KP,1)
        u_KP_history = [u_KP_history nCells * u_UA_KP(transpose(y_KP(i,:)))];
    end
    figure(u_pole_history);
    plot(tSimulation,u_KP_history); hold on;
end
figure(T_pole_history);
ylabel('$T_1$ (K)');
xlabel('Tempo (s)');
legend('LQR', 'Real LQR', '0.9 LQR dominante', '1.1 LQR dominante', '1.3 LQR dominante')
%legend('0.7', '0.9', 'LQR', '1.1', '1.3')
saveas(gcf,'poles_T.eps','epsc')


figure(u_pole_history);
ylabel('$U_A$ (W/K)');
xlabel('Tempo (s)');
legend('LQR', 'Real LQR', '0.9 LQR dominante', '1.1 LQR dominante', '1.3 LQR dominante')
%legend('0.7', '0.9', 'LQR', '1.1', '1.3')
saveas(gcf,'poles_U.eps','epsc')


%%%% Observador %%%%
%%

% @M, here we need realistic values, we can discuss
Vd = 0.1*eye([6,6]);
Vn = 0.7*eye([nMeasuresC,nMeasuresC]);

rng(1,"twister");
disturbance_noise = sqrt(Vd)*randn(6,size(tSimulation,2));
% É equivalente ao ruído gaussiano no intervalo +/- 3C
measure_noise = sqrt(Vn)*randn(2,size(tSimulation,2));

% Criar observador, LQE, com LQR regulador
sys_observer = ss(A_cl_lqr, B, C, D);
L_lqr_obs = (lqr(A_cl_lqr',C',Vd,Vn))';
eig(A_cl_lqr-L_lqr_obs*C);
%%
T_current = TarrSimulation;
T_current_lqr = T_current;
% Use a ordem de grandeza de curso estável como estimativa inicial para aproximação do observador
% A estimativa inicial afeta fortemente a convergência para o valor correto
T_approx_prev = T_equilibrium; 

Tsol = [];
Tsol_approx = [];
Tlqr = [];
Tmeasures = [];

Tsol_approx = [Tsol_approx T_approx_prev];
Tsol = [Tsol T_current];
Tlqr = [Tlqr T_current_lqr];
Tmeasures = [Tmeasures C*(T_current)];

% Resolvemos a dinâmica não linear completa, mas usamos controle linear
for time = 2:size(tSimulation,2)
    T_measured = C*(T_current)+measure_noise(time);
    T_approx_current = T_approx_prev +  dt * (A*(T_approx_prev-T_equilibrium) + B*(u_UA_LQR(T_approx_prev)-U_A) + L_lqr_obs*((T_measured-C*T_equilibrium) - C*(T_approx_prev-T_equilibrium)));
    T_approx_prev = T_approx_current;

    % Método de Euler para resolver equações diferenciais 
    % (em vez de ode45 usado para a mesma tarefa no restante do código)
    % Isso ocorre porque o Luenberger deve conhecer o estado anterior e
    % queremos adicionar ruído manualmente e ainda simular o sistema não linear.
    T_current = T_current + dt*(disturbance_noise(time) + ...
        dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T_current(1),T_current(2),T_current(3),T_current(4),T_current(5),T_current(6),T_ar,u_UA_LQR(T_approx_current),V_cell,V_motor,cH2O,mdot,rho));
    T_current_lqr = T_current_lqr + dt*(disturbance_noise(time) + ...
        dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T_current_lqr(1),T_current_lqr(2),T_current_lqr(3),T_current_lqr(4),T_current_lqr(5),T_current_lqr(6),T_ar,u_UA_LQR(T_current_lqr),V_cell,V_motor,cH2O,mdot,rho));

    Tsol_approx = [Tsol_approx T_approx_current];
    Tsol = [Tsol T_current];
    Tlqr = [Tlqr T_current_lqr];
    Tmeasures = [Tmeasures T_measured];
end
figure;
plot(tSimulation,Tlqr(4,:));hold on
plot(tSimulation,Tsol_approx(4,:),'--');hold on
plot(tSimulation,Tsol(4,:));hold on
legend("LQR com $C=I$", "LQG reconstruida", "LQG real")
ylabel('$T_4$ (K)');
xlabel('Tempo (s)');
saveas(gcf,'LQG_T4.eps','epsc')

figure;
plot(tSimulation,Tmeasures(1,:),'Color',[0.2 0.2 0.9 0.2]);hold on
plot(tSimulation,Tsol_approx(1,:),'--');hold on
plot(tSimulation,Tsol(1,:));hold on
legend("Medida", "LQG reconstruida", "LQG real")
ylabel('$T_1$ (K)');
xlabel('Tempo (s)');
saveas(gcf,'LQG_T1.eps','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Controle clássico               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ANALISE EM MALHA ABERTA %%%%

% Criar o sistema de espaço de estados com a matriz numérica
sys_classic = ss(A, B, C, D);
% Obtendo a função de transferência em malha aberta
G = tf(sys_classic);

G_formula = C*inv((s*eye(size(A)) - A))*B+D

figure;
step(G)
hold on;
step(G_formula)
legend("G ss-modelo", "G função de transferência")


% Função de transferência G_UA_Tx
G_UA_Tfirst = -G % Função de transferência entre T1 e U_A (Assumindo U_A é a segunda entrada)

% please note, G_UA is multiplied by -1 such that a step response of +1 means
% that we cool down the system by 1C, thus the absolute |T| is plotted


disp('Função de Transferência G_{T1, U_A}:');
disp(G_UA_Tfirst);

figure;
subplot(1,2,1)
rlocus(G_UA_Tfirst);
title("Lugar $G:\Delta U_A \longrightarrow \Delta T_1$",'interpreter','latex')
ylabel("Eixo imaginaria")
xlabel("Eixo real")
subplot(1,2,2)
[T_resp_unscaled, t_temp] = step(G_UA_Tfirst);
T_rescaled = T_equilibrium(1) - T_resp_unscaled
plot(t_temp, T_rescaled)
title("Degrau $G:\Delta U_A \longrightarrow\Delta T_1$",'interpreter','latex')
ylabel("$T_1$ (K)")
xlabel("Tempo (s)")
stepinfo(G_UA_Tfirst)
saveas(gcf,'G_UA_Tfirst.eps','epsc')


%%%% Stability with Routh array %%%%
%%

% compensador de atraso
% Compensator made with arthur
G_comp = 0.14286 *(s+0.1842)/(s+0.02632)
G_compensated = G_comp*G_UA_Tfirst

% Coeficientes do numerador e do denominador
[numerator, denominator] = tfdata(G_compensated, 'v');

figure;
subplot(1,2,1)
rlocus(G_compensated);
title("Lugar $G_C G_{U_A \longrightarrow \Delta T_1}$",'interpreter','latex')
ylabel("Eixo imaginaria")
xlabel("Eixo real")
subplot(1,2,2)
[T_resp_unscaled, t_temp] = step(G_compensated);
T_rescaled = T_equilibrium(1) - T_resp_unscaled
plot(t_temp, T_rescaled)
title("Degrau $G_C G_{U_A \longrightarrow \Delta T_1}$",'interpreter','latex')
ylabel("ΔT_1 (K)")
xlabel("Tempo")
stepinfo(G_compensated)
saveas(gcf,'G_UA_Tfirst.eps','epsc')

%%
% search for maximal gain
syms Kmax;
combinedTFpoly = Kmax*numerator + denominator;

% Construindo a tabela de Routh
combinedTFpoly = combinedTFpoly(find(combinedTFpoly, 1):end); % Remove leading zeros if any
n = length(combinedTFpoly);
m = ceil(n / 2);
Routh = Kmax*zeros(n, m);

Routh(1, :) = combinedTFpoly(1:2:end); % Coeficientes de termos ímpares
Routh(2, 1:length(combinedTFpoly(2:2:end))) = combinedTFpoly(2:2:end); % Coeficientes de termos pares

for i = 3:n
    for j = 1:(m-1)
        Routh(i, j) = (Routh(i-1, 1) * Routh(i-2, j+1) - Routh(i-2, 1) * Routh(i-1, j+1)) / Routh(i-1, 1);
    end
    if all(Routh(i, :) == 0) % Special case: row of zeros
        Routh(i, :) = (n-i) * polyder(Routh(i-1, :));
    end
    if Routh(i, 1) == 0 % Avoid division by zero
        Routh(i, 1) = 1e-6; % Small value instead of zero
    end
end

% Solve for ultimate gain by increasing the value of K

for i = 3:n
    % Use a large initial value to get positive solutions. Requieres manual
    % tweaking
    fzero(matlabFunction(Routh(i,1)),1000)
end
% Manual inspection of intersections about to find ultimate gain that was
% printed above. Se on Routh array if values actually is ultimate
Kmax = 929.9327
disp("Com Kultimate")
disp(double(subs(Routh)));
Kmax = Kmax * 1.01
disp("Com 101% Kultimate")
disp(double(subs(Routh)));
Kmax = Kmax * 0.99 / 1.01;
disp("Com 99% Kultimate")
disp(double(subs(Routh)));
Kmax = Kmax/0.99;

Kultimate = Kmax;

% Analisando a estabilidade malha aberta
Kmax = 0
RouthAberta = subs(Routh)

disp("Stability binary, 1) According to Routh calculation 2) Matlab $ stable")
stability = all(RouthAberta(:, 1) > 0)
disp(isstable(G_compensated))

disp('Tabela de Routh Malha Aberta:');
disp(double(RouthAberta));

%%
RouthUltimate = subs(Routh);
stability = all(RouthUltimate(:, 1) > 0)
TF_first_ultimate = feedback(Kultimate * G_compensated,1)
disp(isstable(TF_first_ultimate))

disp('Tabela de Routh Malha Kultimate:');
disp(double(RouthUltimate));

cl_poles_u = pole(TF_first_ultimate)
osc_arr = imag(sort(cl_poles_u,'ComparisonMethod','real'))
P_u = 2*pi / abs(osc_arr(end))

% Calculate PID parameters
Kp_PID = 0.6 * Kultimate;
Ti_PID = 0.5 * P_u;
Td_PID = 0.125 * P_u;

% PI controller (proportional and integral action)
Kp_PD = 0.8 * Kultimate;
Ti_PD = 0.125 * P_u;

% P controller (no integral or derivative action)
Kp_P = 0.5 * Kultimate;

PID_controller = pid(Kp_PID,Kp_PID/Ti_PID,Kp_PID*Td_PID)
PD_controller = pid(Kp_PD,0, Kp_PD*Ti_PD)
P_controller = pid(Kp_P,0,0)

PID_sys = feedback(PID_controller* G_compensated,1)
PD_sys = feedback(PD_controller* G_compensated,1)
P_sys = feedback(P_controller* G_compensated,1)

% PI-controller manually tuned
G_PI_manual = 5.2208 * (1+4*s)/s 

PI_manual_sys = feedback(G_PI_manual* G_UA_Tfirst,1)

G_lead_comp = 55.892*(s+1.476)/(s+5.056)
lead_comp_sys = G_lead_comp*G_UA_Tfirst

step_time = 60; % s

figure;
[T_resp_unscaled, t_temp] = step(PID_sys,step_time );
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

[T_resp_unscaled, t_temp] = step(PD_sys,step_time);
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

[T_resp_unscaled, t_temp] = step(P_sys,step_time);
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

[T_resp_unscaled, t_temp] = step(G_compensated,step_time);
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

title("Degrau $G_C G_{U_A \longrightarrow \Delta T_1}$",'interpreter','latex')
ylabel("$T_1$ (K)")
xlabel("Tempo (s)")
legend("PID, ZN", "PD, ZN", "P, ZN","So compensador")
saveas(gcf,'G_comp_ZN_PID.eps','epsc')

figure;
G_ZN_PID_U = PID_sys/(1+PID_sys*G_compensated)
[U_A_resp_unscaled,t_resp] = step(G_ZN_PID_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

G_ZN_PD_U = PD_sys/(1+PD_sys*G_compensated)
[U_A_resp_unscaled,t_resp] = step(G_ZN_PD_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

G_ZN_P_U = P_sys/(1+P_sys*G_compensated)
[U_A_resp_unscaled,t_resp] = step(G_ZN_P_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

G_ZN_sem_U = 1/(1+1*G_compensated)
[U_A_resp_unscaled,t_resp] = step(G_ZN_sem_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

legend("PID, ZN", "PD, ZN", "P, ZN","So compensador")

ylabel("$U_A$ (W/K)")
xlabel("Tempo (s)")

figure;

[T_resp_unscaled, t_temp] = step(PI_manual_sys ,step_time );
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

[T_resp_unscaled, t_temp] = step(feedback(lead_comp_sys,1),step_time);
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

[T_resp_unscaled, t_temp] = step(G_UA_Tfirst,step_time);
T_rescaled = T_equilibrium(1) - T_resp_unscaled;
plot(t_temp, T_rescaled);hold on

title("Degrau $G_{U_A \longrightarrow \Delta T_1}$",'interpreter','latex')
ylabel("$T_1$ (K)")
xlabel("Tempo (s)")
legend("Ajustado manualmente", "lead with closed", "aberta")
saveas(gcf,'G_comp_ZN_PID_U.eps','epsc')


figure;
G_manual_PID_U = PI_manual_sys/(1+PI_manual_sys*G_UA_Tfirst)
[U_A_resp_unscaled,t_resp] = step(G_manual_PID_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

G_lead_U = lead_comp_sys/(1+lead_comp_sys*G_UA_Tfirst)
[U_A_resp_unscaled,t_resp] = step(G_lead_U,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

G_aberta = G_UA_Tfirst
[U_A_resp_unscaled,t_resp] = step(G_aberta,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;
plot(t_resp,U_A_resp_scaled); hold on

legend("Ajustado manualmente", "lead with closed", "aberta")

ylabel("$U_A$ (W/K)")
xlabel("Tempo (s)")



% Since P_u = inf, we cannot use standard Z-N-tuning
% indeed looking at the step response we have τdead = 0
% which means Ziegler-Nichols Open-Loop Tuning Method / Process Reaction Method
% also will not apply. We can used P-controller with obtained values
% however although it is expected to be bad

%%%% Malha Fechada %%%%
%%
figure;
step(G_UA_Tfirst)
hold on;
step(feedback(G_UA_Tfirst,1))
%legend("Malha aberta", "Malha fechada")

%%

s = tf("s")




% H = disturbance


G_UA_Tfirst_C = feedback(G_C*G_UA_Tfirst,1)

figure;
subplot(2,1,1)
[T_resp_unscaled,t_resp] = step(G_UA_Tfirst_C,step_time);
T_resp_scaled = T_equilibrium(1)-T_resp_unscaled;
plot(t_resp,T_resp_scaled)
title("Degrau $(G_CG_{\Delta U_A \longrightarrow \Delta T_1})/(1+G_CG_{\Delta U_A \longrightarrow \Delta T_1})$",'interpreter','latex')
ylabel("$T_1$ (K)")
stepinfo(G_UA_Tfirst_C)
subplot(2,1,2)
G_UA_Tfirst_C_UA = G_C/(1+G_C*G_UA_Tfirst)
[U_A_resp_unscaled,t_resp] = step(G_UA_Tfirst_C_UA,step_time); % plot actual value during response
U_A_resp_scaled = U_A*nCells + nCells*U_A_resp_unscaled;

disp("Final value")
disp(U_A_resp_scaled(end))
disp("Max value")
disp(max(U_A_resp_scaled))
disp("Min value")
disp(min(U_A_resp_scaled))

plot(t_resp,U_A_resp_scaled)
resp_info = stepinfo(t_resp, U_A_resp_scaled)
ylabel("$U_A$ (W/K)")
xlabel("Tempo (s)")
saveas(gcf,'G_UA_Tfirst.eps','epsc')
%%
figure;
%subplot(2,2,1)
rlocus(G_UA_Tfirst_C);
title("Lugar $G_C*G_{\Delta U_A \longrightarrow \Delta T_1}/(1+G_C*G_{\Delta U_A \longrightarrow \Delta T_1})$",'interpreter','latex')
ylabel("Eixo imaginaria")
xlabel("Eixo real")
%%


% P ZN - tuning
%K_c_zn = 0.6 * Kultimate
%tau_i = P_u/2
%tau_d = P_u/8
%G_c_zn_PID = K_c_zn * (1 + tf([1], [tau_i 0]) + tf([tau_d 0], [tau_d/N_d 1]))

%%
pidCoef_ZN_PID = 1.4%pid(0.5 * Kultimate,0.6 * Kultimate,0.01 * Kultimate)
pidCoef_ZN_PI = 1.2%pid(1.5 * Kultimate,0.6 * Kultimate,0.01 * Kultimate)
pidCoef_ZN_P = 1%pid(2.5 * Kultimate,0.6 * Kultimate,0.01 * Kultimate)

G_UA_Tfirst_PID_ZN = feedback(pidCoef_ZN_PID*G_UA_Tfirst,1)
G_UA_Tfirst_PI_ZN = feedback(pidCoef_ZN_PI*G_UA_Tfirst,1)
G_UA_Tfirst_P_ZN = feedback(pidCoef_ZN_P*G_UA_Tfirst,1)


pidCoef_manual = 0.2%pid(0.5 * Kultimate,0.01 * Kultimate,0.01 * Kultimate)

G_UA_Tfirst_PID_manual = feedback(pidCoef_manual*G_UA_Tfirst,1)


%%
figure;
step(G_UA_Tfirst_PID_manual); hold on;
step(G_UA_Tfirst_PID_ZN); hold on;
step(G_UA_Tfirst_PI_ZN); hold on;
step(G_UA_Tfirst_P_ZN); hold on;
step(G_UA_Tfirst);
ylabel("T_1 (K)")
xlabel("Tempo (s)")
legend("Manual PID","ZN PID", "ZN PI", "ZN P", "Malha aberta")
saveas(gcf,'PID_step_response.eps','epsc')

%%

figure;
subplot(1,2,1)
rlocus(G_UA_Tend);
title("Lugar $G:\Delta U_A \longrightarrow \Delta T_6$",'interpreter','latex')
ylabel("Eixo imaginaria")
xlabel("Eixo real")
subplot(1,2,2)
step(G_UA_Tend);
title("Degrau $G:\Delta U_A \longrightarrow \Delta T_6$",'interpreter','latex')
ylabel("ΔT_6 (K)")
xlabel("Tempo")
stepinfo(G_UA_Tend)
saveas(gcf,'G_UA_Tend.eps','epsc')


% TODO: implement H

%% Code de Arthur %%%


% Função de transferência G_UA
% G_UA_original = G(1, 1); % Função de transferência entre T1 e U_A 

% Coeficientes do numerador e do denominador da função de transferência original
% [numerator, denominator] = tfdata(G_UA_original, 'v');

% Multiplicar a função de transferência por -1 para ajuste correto
G_UA = G_UA_Tfirst

% Mostrar a função de transferência ajustada
disp('Função de Transferência G_{T1, U_A} * -1:');
disp(G_UA);


% Encontre os polos e zeros
poles = pole(G_UA);
zeros = zero(G_UA);

% Exiba os valores dos polos e zeros
disp('Polos:');
disp(poles);
disp('Zeros:');
disp(zeros);

% Plotar o Root Locus
figure;
rlocus(G_UA);
title('Root Locus de G_{T1, U_A}');
grid on;

% Definição da faixa de frequência desejada
fmin = 10^-2;  % Frequência mínima
fmax = 10^1;   % Frequência máxima
w = logspace(log10(fmin), log10(fmax), 500);  % Vetor de frequência

% Plotar o diagrama de Bode na faixa de frequência especificada
figure;
bode(G_UA, w);
grid on;
title('Diagrama de Bode de G_UA com faixa de frequência ajustada');

% Analisar a resposta ao degrau
figure;
step(G_UA);
title('Resposta ao Degrau de G_{T1, U_A}');


%%% CONTROLADOR P, PI, PD e PID %%%

% Parâmetros iniciais do controlador
Ti = 999; % Definir o tempo integral para um valor muito alto
Td = 0; % Definir o tempo derivativo para zero

% Ajustar ganho proporcional Kc
Ku = 0; % Inicialização do ganho final
Pu = 0; % Inicialização do período de oscilação
Kc_values = 0:0.1:30; % Valores de Kc a serem testados

% Loop para encontrar Ku e Pu
found = false;
for Kc = Kc_values
    PID_Controller = pid(Kc, Kc/Ti, Kc*Td);
    sys_cl = feedback(PID_Controller * G_UA, 1);
    t = 0:0.01:100;
    y = step(sys_cl, t);
    
    % Verificar oscilações
    [pks, locs] = findpeaks(y); % Encontrar picos na resposta
    
    if length(pks) > 2
        % Calcular o período de oscilação
        Pu = mean(diff(t(locs)));
        Ku = Kc;
        found = true;
        break;
    end
end

% Verificar se Ku e Pu foram encontrados
if ~found
    error('Não foi possível encontrar Ku e Pu. Ajuste o intervalo de Kc e tente novamente.');
else
    disp(['Ku encontrado: ', num2str(Ku)]);
    disp(['Pu encontrado: ', num2str(Pu)]);
end

% Parâmetros PID usando Ziegler-Nichols
Kc = Ku / 1.7;
Ti = Pu / 2;
Td = Pu / 8;

% Exibir os valores calculados
disp('Parâmetros PID calculados:')
disp(['Kc: ', num2str(Kc)])
disp(['Ti: ', num2str(Ti)])
disp(['Td: ', num2str(Td)])

% Criar controladores P, PI, PD e PID
P_Controller = pid(Kc, 0, 0);
PI_Controller = pid(Kc, Kc/Ti, 0);
PD_Controller = pid(Kc, 0, Kc*Td);
PID_Controller = pid(Kc, Kc/Ti, Kc*Td);

% Definir os sistemas em malha fechada para cada controlador
sys_cl_P = feedback(P_Controller * G_UA, 1);
sys_cl_PI = feedback(PI_Controller * G_UA, 1);
sys_cl_PD = feedback(PD_Controller * G_UA, 1);
sys_cl_PID = feedback(PID_Controller * G_UA, 1);

% Simulação de resposta ao degrau
t = 0:0.01:60;
[y_P, t] = step(sys_cl_P, t);
[y_PI, t] = step(sys_cl_PI, t);
[y_PD, t] = step(sys_cl_PD, t);
[y_PID, t] = step(sys_cl_PID, t);

% Plotar a resposta ao degrau para cada controlador
figure;
plot(t, y_P, 'r', 'DisplayName', 'Controlador P');
hold on;
plot(t, y_PI, 'g', 'DisplayName', 'Controlador PI');
plot(t, y_PD, 'b', 'DisplayName', 'Controlador PD');
plot(t, y_PID, 'k', 'DisplayName', 'Controlador PID');
hold off;
title('Resposta ao Degrau com Diferentes Tipos de Controladores')
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('show');
grid on;

% Avaliar o overshoot e o tempo de estabilização para cada controlador
info_P = stepinfo(y_P, t);
overshoot_P = info_P.Overshoot;
settling_time_P = info_P.SettlingTime;
disp(['Overshoot Controlador P: ', num2str(overshoot_P), ' %']);
disp(['Tempo de Estabilização Controlador P: ', num2str(settling_time_P), ' s']);

info_PI = stepinfo(y_PI, t);
overshoot_PI = info_PI.Overshoot;
settling_time_PI = info_PI.SettlingTime;
disp(['Overshoot Controlador PI: ', num2str(overshoot_PI), ' %']);
disp(['Tempo de Estabilização Controlador PI: ', num2str(settling_time_PI), ' s']);

info_PD = stepinfo(y_PD, t);
overshoot_PD = info_PD.Overshoot;
settling_time_PD = info_PD.SettlingTime;
disp(['Overshoot Controlador PD: ', num2str(overshoot_PD), ' %']);
disp(['Tempo de Estabilização Controlador PD: ', num2str(settling_time_PD), ' s']);

info_PID = stepinfo(y_PID, t);
overshoot_PID = info_PID.Overshoot;
settling_time_PID = info_PID.SettlingTime;
disp(['Overshoot Controlador PID: ', num2str(overshoot_PID), ' %']);
disp(['Tempo de Estabilização Controlador PID: ', num2str(settling_time_PID), ' s']);

% Plotar o Root Locus
figure;
rlocus(sys_cl_P);
title('Root Locus de G_{T1, U_A} com P');
grid on;

figure;
rlocus(sys_cl_PI);
title('Root Locus de G_{T1, U_A} com PI');
grid on;

figure;
rlocus(sys_cl_PD);
title('Root Locus de G_{T1, U_A} com PD');
grid on;

figure;
rlocus(sys_cl_PID);
title('Root Locus de G_{T1, U_A} com PID');
grid on;

% Encontre os polos e zeros para análise apenas do PID
poles1 = pole(sys_cl_PID);
zeros1 = zero(sys_cl_PID);

% Exiba os valores dos polos e zeros
disp('Polos PID:');
disp(poles1);
disp('Zeros PID:');
disp(zeros1);

% Plotar diagramas de Bode
figure;
bode(G_UA, w);
hold on;
bode(sys_cl_P, w);
title('Diagrama de Bode: Controlador P vs Malha Aberta');
legend('Malha Aberta', 'Controlador P');
grid on;

figure;
bode(G_UA, w);
hold on;
bode(sys_cl_PI, w);
title('Diagrama de Bode: Controlador PI vs Malha Aberta');
legend('Malha Aberta', 'Controlador PI');
grid on;

figure;
bode(G_UA, w);
hold on;
bode(sys_cl_PD, w);
title('Diagrama de Bode: Controlador PD vs Malha Aberta');
legend('Malha Aberta', 'Controlador PD');
grid on;

figure;
bode(G_UA, w);
hold on;
bode(sys_cl_PID, w);
title('Diagrama de Bode: Controlador PID vs Malha Aberta');
legend('Malha Aberta', 'Controlador PID');
grid on;

%%% MARGEM DE GANHO E MARGEM DE FASE %%%

figure;
margin(G_UA, w);
title('Margem de Ganho e Margem de Fase em malha aberta');
grid on;

[GM, PM, Wcg, Wcp] = margin(G_UA);
disp(['Margem de Ganho: ', num2str(GM), ' dB']);  % Margem infinita
disp(['Margem de Fase: ', num2str(PM), ' graus']);  % Nunca chega aos 180 graus
disp(['Frequência de Cruzamento de Ganho: ', num2str(Wcg), ' rad/s']);
disp(['Frequência de Cruzamento de Fase: ', num2str(Wcp), ' rad/s']);


%%% COMPENSADOR DE AVANCO ATRAVES DO DIAGRAMA DE BODE %%%

% Determinar a Margem de Fase Desejada
desired_phase_margin = 45; % Margem de fase desejada em graus, NAO SEI QUAL O VALOR CERTO
additional_phase_margin = desired_phase_margin - PM;

% Projetar o Compensador de Avanço
% Determinar a frequência onde o ganho cruzado acontece
freq = Wcp; % Utilizando a frequência de cruzamento de fase, TALVEZ O PROBLEMA SEJA AQUI, FREQUENCIA BAIXA

% Calcular o ângulo necessário de avanço (em radianos)
phi = (additional_phase_margin) * (pi/180);

% Determinar alpha e T
alpha = (1 - sin(phi)) / (1 + sin(phi));
T = 1 / (freq * sqrt(alpha));  % Usando Wcp para determinar T

% Criação do compensador de avanço
C = tf([T 1], [alpha*T 1]);

% Implementar o Compensador
sys_comp = series(C, G_UA);

% Nova Resposta em Frequência
figure;
bode(G_UA, sys_comp, w);
legend('Malha aberta', 'Sistema com Compensador de Avanço', 'Location', 'southeast');
grid on;
title('Diagrama de Bode com e sem Compensador de Avanço');

% Análise de Margem de Fase
figure;
margin(sys_comp, w);
title('Margem de Fase do Sistema com Compensador de Avanço');

% Resposta ao degrau
figure;
step(sys_comp, 0:0.01:1200);
legend('Sistema com Compensador de Avanço', 'Location', 'southeast');
title('Resposta ao Degrau com Compensador de Avanço');
grid on;

% Avaliar a Resposta ao Degrau
info = stepinfo(sys_comp);
disp(['Overshoot com controlador de avanco: ', num2str(info.Overshoot), ' %']);
disp(['Settling Time com controlador de avanco: ', num2str(info.SettlingTime), ' seconds']);
