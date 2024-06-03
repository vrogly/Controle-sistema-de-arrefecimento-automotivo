clc;clear;
close all
set(groot,'defaultAxesFontSize',16)
set(groot,'DefaultFigureRenderer','painters');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

disp("Running")

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
nCells = nT-1; % We have one more T (surface) than cells


u_params = [U_A]; % Note! You will have do do manual changes if you increase params here by adding more u_UA etc functions
fluctating_params = [Qdotger, mdot,Delta_P_H2O,T_ar];

A = jacobian(dTi, Tarr);
B = jacobian(dTi, u_params);
E = jacobian(dTi, fluctating_params);

latexify("\frac{d}{dt}\vec T =",dTi)
latexify("A=",A)
latexify("B=",B)

% y independent of x after a substitution, i.e. we cant update the value again

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

F = 0.9;
T_ar = 298; % K
U_A = 100/nCells; % W/K

dTi = subs(dTi);

% Simultion used for comparision
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

nMeasuresC = 2;
C = zeros([nMeasuresC, nT]); % We only measure 2 values
C(1) = 1; % We can only measure edge temperatures
C(end) = C(1);
for nn = 2:nMeasuresC-1
    C(nn,nn) = 1;
end

D_full_observe = zeros([size(C_full_observe,1),size(B,2)]); % Matriz de zeros
D = zeros([size(C,1),size(B,2)]); % Matriz de zeros

% Definindo o sistema
sys_ol = ss(A, B, C_full_observe, D_full_observe)

% Obtendo os polos
polos = eig(A);
disp('Polos:');
disp(polos);

% Plotando o diagrama de Bode
figure;
bode(sys_ol);
title('Diagrama de Bode');


%%%% REGULADOR LINEAR QUADRÁTICO %%%%

% Definindo as matrizes de peso Q e R
Q = diag([1,0,0,0,0,1]);%
% Did not help to put 0.1 or 0.5 instead of zeros, only made worse 
R = 0.01; % Peso para a entrada, %% Make it expensive to change fan speed

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

%% Compare open loop to closed loop with LQR

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





%% Pole Placement

% Uses LQR poles and tries to vary them slightly to manually
% decide if there are more suitable placements
% plot the control u over time to see that it does not exceed any
% specifications

polesFactorArr = [0.7, 0.9, 1, 1.1, 1.3];
T_pole_history = figure; 
u_pole_history = figure;

for poleFactor = 1:size(polesFactorArr,2)
    desired_poles = polos_lqr * polesFactorArr(poleFactor)
    Kp = place(A, B, desired_poles);
    A_poleTemp = A - B * Kp;
    disp(eig(A_poleTemp))
    
    % Here we could build and simulate a state space model like
    % sys_poleTemp = ss(A_poleTemp, B, C_full_observe, D_full_observe);
    % But instead we simulate full non-linear dynamics and update the
    % control vairable using our linear model

    u_UA_KP = @(T)double(U_A-Kp*(T-T_equilibrium));
    [t,y_KP] = ode45(@(t,T)dT_pure(Delta_P_R,Delta_P_H2O,F,Qdotger,T(1),T(2),T(3),T(4),T(5),T(6),T_ar,u_UA_KP(T),V_cell,V_motor,cH2O,mdot,rho),tSimulation,TarrSimulation);    
    
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
%title("Controlling initial conditions to desired state")
ylabel('$T_1$ (K)');
xlabel('Tempo (s)');
legend('0.7', '0.9', 'LQR', '1.1', '1.3')
saveas(gcf,'poles_T.eps','epsc')


figure(u_pole_history);
%title("Controlling initial conditions to desired state")
ylabel('$U_A$ (W/K)');
xlabel('Tempo (s)');
legend('0.7', '0.9', 'LQR', '1.1', '1.3')
saveas(gcf,'poles_U.eps','epsc')


%% Luenberger

% @M, here we need realistic values, we can discuss
Vd = 0.1*eye([6,6]);
Vn = 0.1*eye([2,2]);

rng(1,"twister");
disturbance_noise = sqrt(Vd)*randn(6,size(tSimulation,2));
measure_noise = sqrt(Vn)*randn(2,size(tSimulation,2));


% Create observer, LQE, with LQR regulator
sys_observer = ss(A_cl_lqr, B, C, D);
L_lqr_obs = (lqr(A_cl_lqr',C',Vd,Vn))';
eig(A_cl_lqr-L_lqr_obs*C);

T_current = TarrSimulation;
T_current_lqr = T_current;
% Use the order of magnitude as an initial guess,
% the initial guess heavily effects the convergence
% to the correct value
T_approx_prev = T_equilibrium; 

Tsol = [];
Tsol_approx = [];
Tlqr = [];
Tmeasures = [];

Tsol_approx = [Tsol_approx T_approx_prev];
Tsol = [Tsol T_current];
Tlqr = [Tlqr T_current_lqr];
Tmeasures = [Tmeasures C*(T_current)];

% We solve for full non linear dynamics but use linear control
for time = 2:size(tSimulation,2)
    T_measured = C*(T_current)+measure_noise(time);
    T_approx_current = T_approx_prev +  dt * (A*(T_approx_prev-T_equilibrium) + B*(u_UA_LQR(T_approx_prev)-U_A) + L_lqr_obs*((T_measured-C*T_equilibrium) - C*(T_approx_prev-T_equilibrium)));
    T_approx_prev = T_approx_current;

    % Euler's method to solve differential equations (instead of ode45 as
    % used for the same task in the rest of the code. This is beacuse the
    % Luenberger must know previous state and we want to manually add noise
    % and still simulate the non linear system.
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

%%

%%% ANALISE EM MALHA ABERTA %%%

% Criar o sistema de espaço de estados com a matriz numérica
sys = ss(A, B, C, D);

% Obtendo a função de transferência em malha aberta
G = tf(sys);

% Função de transferência G_UA
G_UA = G(1, 1); % Função de transferência entre T1 e U_A (Assumindo U_A é a segunda entrada)

% Coeficientes do numerador e do denominador
[numerator, denominator] = tfdata(G_UA, 'v');

% Criando a função de transferência
G_UA = tf(numerator, denominator);
disp('Função de Transferência G_{T1, U_A}:');
disp(G_UA);

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

inEqArr = [];
for i = 3:n
    if isempty(solve((Routh(i,1)) > 0)) == 0
        inEqArr = [inEqArr ((Routh(i,1)) > 0)];
    else
        disp(simplify(Routh(i,1)))
    end
end

inEqArr

solve(inEqArr,[Kmax],'ReturnConditions',true)

% Manually insert the obtained value
Kultimate = 98862924258641/329469785125776


% Analisando a estabilidade malha aberta
Kmax = 0
RouthAberta = subs(Routh)
stability = all(RouthAberta(:, 1) > 0)
G_UA
disp(isstable(G_UA))

disp('Tabela de Routh Malha Aberta:');
disp(double(RouthAberta));


Kmax = Kultimate;
RouthUltimate = subs(Routh);
stability = all(RouthUltimate(:, 1) > 0)
TF_ultimate = feedback(G_UA,Kmax)
disp(isstable(TF_ultimate))

disp('Tabela de Routh Malha Kultimate:');
disp(double(RouthUltimate));

cl_poles_u = pole(TF_ultimate)
P_u = 2*pi / imag(cl_poles_u(1))
