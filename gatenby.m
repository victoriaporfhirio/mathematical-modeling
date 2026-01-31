clear; clc;

% ================================
% Parâmetros do modelo
% ================================
r1 = 1.1;        % taxa de crescimento das células tumorais (N1)
r2 = 2.2;        % taxa de crescimento das células normais (N2)
K1 = 10;         % capacidade de suporte das células tumorais
K2 = 10;         % capacidade de suporte das células normais
alpha12 = 0.005; % efeito de N2 sobre N1
alpha21 = 0.1;   % efeito de N1 sobre N2

% ================================
% Condições iniciais
% ================================
t0 = 0.0;        % tempo inicial (dias)
tf = 40;         % tempo final (dias)
N1(1) = 5;       % população inicial das células tumorais
N2(1) = 8;       % população inicial das células normais

% ================================
% Configuração do método de Euler
% ================================
h = 0.01;                          % passo de integração (variações testadas: 0.5, 0.1, 0.01, 0.001)
n_iter = round((tf - t0)/h);    % número de iterações

% ================================
% Iterações do método de Euler
% ================================
for n = 1:n_iter
    N1(n+1) = N1(n) + h * r1 * N1(n) * (1 - (N1(n) + alpha12 * N2(n)) / K1);
    N2(n+1) = N2(n) + h * r2 * N2(n) * (1 - (N2(n) + alpha21 * N1(n)) / K2);
end

% ================================
% Vetor de tempo
% ================================
t = t0:h:tf;

% ================================
% Gráfico das populações
% ================================
figure;
plot(t, N1, 'r-', 'LineWidth', 2); hold on;
plot(t, N2, 'b-', 'LineWidth', 2);
xlabel('Tempo (dias)');
ylabel('População');
legend('Células Tumorais (N1)', 'Células Normais (N2)');
title('Modelo de Gatenby - Competição entre células tumorais e normais');
grid on;

