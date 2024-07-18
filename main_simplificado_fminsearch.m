% Dados experimentais
t_exp = [2.5, 5.0, 10.0, 15.0, 30.0, 45.0, 60.0, 75.0]; % em minutos
C_H2_exp = [1.4, 1.4, 1.8, 1.8, 2.4, 2.0, 2.6, 3.2]; % Concentrações experimentais (em mmol/g)
C_CO_exp = [2.8, 3.0, 2.8, 3.0, 1.2, 1.5, 1.95, 1.6];
C_CH4_exp = [4.1, 4.8, 5.1, 5.2, 5.6, 5.0, 5.9, 6.3];
C_CO2_exp = [4.4, 4.7, 5.7, 5.2, 6.2, 5.1, 6.0, 7.6];
C_exp = [C_H2_exp; C_CO_exp; C_CH4_exp; C_CO2_exp]'; % Matriz de concentrações experimentais
C0 = [5.095, 0, 0, 0, 0, 0, 0]; % Condições iniciais
tspan = [0 75]; % Intervalo de tempo em minutos
T = 600;

x = 10; y = 12; z = 4; % Definição de x, y, z

% Chute inicial para as constantes de reação
k = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];

%% Otimização
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1000);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[k_opt, error_opt] = fminsearch(@(k) objective_function(k, tspan, C0, t_exp, x, y, z, C_exp), k, options);

% Recalcular as concentrações com as constantes otimizadas
[t_opt, C_opt] = ode45(@(t, C) reactions(t, C, k_opt, x, y, z), tspan, C0);

%% Plotar os resultados
figure;
hold on;
plot(t_exp, C_exp(:, 1), 'ro', 'DisplayName', 'H2 Experimental');
plot(t_exp, C_exp(:, 2), 'go', 'DisplayName', 'CO Experimental');
plot(t_exp, C_exp(:, 3), 'bo', 'DisplayName', 'CH4 Experimental');
plot(t_exp, C_exp(:, 4), 'ko', 'DisplayName', 'CO2 Experimental');
plot(t_opt, C_opt(:, 2), 'r-', 'DisplayName', 'H2 Simulado');
plot(t_opt, C_opt(:, 3), 'g-', 'DisplayName', 'CO Simulado');
plot(t_opt, C_opt(:, 4), 'b-', 'DisplayName', 'CH4 Simulado');
plot(t_opt, C_opt(:, 5), 'k-', 'DisplayName', 'CO2 Simulado');
xlabel('Tempo (min)');
ylabel('Concentração (mmol/g)');
legend show;
hold off;

%% Função para resolver as ODEs
function dCdt = reactions(t, C, k, x, y, z)
    % Constantes para reações
    K10 = 2.68;
    K11 = 1.02e4;
    k(12) = k(10)/K10;
    k(13) = k(11)/K11;
    
    % Equações diferenciais
    dCdt = zeros(7, 1);
    dCdt(1) = -k(2)*C(1);
    dCdt(2) = (x-z+y/2)*k(3)*C(7)*C(6) + (2*x-z+y/2)*k(4)*C(7)*C(6) + k(8)*C(7) + k(10)*C(3)*C(6)...
        - k(12)*C(4)*C(2) - 3*k(11)*C(3)*C(2) + 3*k(13)*C(5)*C(6);
    dCdt(3) = x*k(3)*C(7)*C(6) + k(5)*C(7) - k(10)*C(3)*C(6) + k(12)*C(4)*C(2) - k(11)*C(3)*C(2)...
        + k(13)*C(5)*C(6);
    dCdt(4) = x*k(4)*C(7)*C(6) + k(6)*C(7) + k(10)*C(3)*C(6) - k(12)*C(4)*C(2);
    dCdt(5) = k(7)*C(7) + k(11)*C(3)*C(2) - k(13)*C(5)*C(6);
    dCdt(6) = -(2*x-z)*k(4)*C(7)*C(6) - k(10)*C(3)*C(6) + k(12)*C(4)*C(2) + k(11)*C(3)*C(2) - k(13)*C(5)*C(6);
    dCdt(7) = k(2)*C(1) - k(3)*C(7)*C(6) - k(4)*C(7)*C(6) - k(9)*C(7);
end

%% Função objetivo para otimização
function error = objective_function(k, tspan, C0, t_exp, x, y, z, C_exp)
    % Resolver as ODEs com as constantes k
    [t, C] = ode45(@(t, C) reactions(t, C, k, x, y, z), tspan, C0);
    
    % Interpolar os resultados simulados nos tempos experimentais
    C_sim_interp = zeros(length(t_exp), size(C_exp, 2));
    for i = 1:size(C_exp, 2)
        C_sim_interp(:, i) = interp1(t, C(:, i+1), t_exp); % Começando de i+1 para ignorar a lignina
    end
    
    % Calcular o erro quadrático
    error = sum(sum((C_sim_interp - C_exp).^2)); % Ignorar a lignina na comparação
end