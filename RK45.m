clear; clc; close all;

% Parâmetros iniciais
y_an = 100;  % População inicial (analítica)
y_45 = 100;  % População inicial (RK45)
t = 0;       % Tempo inicial
T = 30;      % Tempo final
r = 0.15;    % Taxa de crescimento populacional

% Parâmetros do controle adaptativo
hMax = 1;          % Passo máximo
hMin = 1e-7;       % Passo mínimo
h = hMax;          % Passo inicial
tol = 1e-2;        % Tolerância do erro

% Vetores para armazenar resultados
t_values = t;
y_an_values = y_an;
y_45_values = y_45;

while true
    N = y_45_values(end); % População atual

    % Cálculo dos coeficientes do RK45 (Runge–Kutta–Fehlberg)
    k1 = h * r * N;
    k2 = h * r * (N + k1/4);
    k3 = h * r * (N + (3/32)*k1 + (9/32)*k2);
    k4 = h * r * (N + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    k5 = h * r * (N + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    k6 = h * r * (N - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    % Estimativa do erro
    R = (1/h) * abs((1/360)*k1 - (128/4275)*k3 - (2197/75240)*k4 + (1/50)*k5 + (2/55)*k6);

    % Se o erro for aceitável, avança o passo
    if R <= tol
        t_next = t_values(end) + h;
        t_values(end+1) = t_next;

        % Solução analítica
        y_an_values(end+1) = y_an_values(1) * exp(r * t_next);

        % Atualiza RK45
        r45 = (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;
        y_45_values(end+1) = y_45_values(end) + r45;
    end

    % Controle adaptativo do passo
    q = 0.84 * (tol / R)^(0.25);

    if q <= 0.1
        h = h / 10;
    elseif q >= 4
        h = 4 * h;
    else
        h = q * h;
    end

    % Limites do passo
    if h > hMax
        h = hMax;
    end

    % Critérios de parada
    if t_values(end) >= T
        break;
    elseif t_values(end) + h > T
        h = T - t_values(end);
    elseif h < hMin
        warning('Passo mínimo atingido. Interrompendo simulação.');
        break;
    end
end

% Plotagem dos resultados
figure;
plot(t_values, y_an_values, 'k-', 'LineWidth', 1.2); hold on;
plot(t_values, y_45_values, 'r--', 'LineWidth', 1.2);
xlabel('Instante de tempo');
ylabel('População');
legend('Solução Analítica', 'Método RK45');
title('Crescimento Populacional - Método de Runge–Kutta–Fehlberg (RK45)');
grid on;

% Cálculo da norma da diferença entre as soluções
dif_norma = norm(y_an_values - y_45_values);

% Exibição dos resultados
fprintf('Norma da diferença entre analítica e RK45: %.10f\n', dif_norma);
fprintf('Número de passos realizados: %d\n', length(y_an_values)-1);

