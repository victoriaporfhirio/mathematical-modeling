% Método de Euler e Diferenças Concentradas para o crescimento exponencial

% Parâmetros do problema
y0 = 10;       % condição inicial
r = 0.15;      % taxa de crescimento
t0 = 0;        % tempo inicial
tf = 10;       % tempo final

% Passos a serem testados
hs = [1, 0.5, 0.1];

% Figura para o gráfico
figure;
hold on;

% Guardar "handles" das curvas para legenda
handles = [];
labels = {};

% Loop para cada h
for j = 1:length(hs)
    h = hs(j);
    N = (tf - t0)/h;

    % Vetores de tempo
    t = t0:h:tf;

    % ---------------------------
    % Método de Euler
    % ---------------------------
    y_euler = zeros(1, N+1);
    y_euler(1) = y0;
    for n = 1:N
        y_euler(n+1) = y_euler(n) + h * r * y_euler(n);
    end
    hplot = plot(t, y_euler, 'o-', 'LineWidth', 1.2);
    handles(end+1) = hplot;
    labels{end+1} = sprintf('Euler (h = %.2f)', h);

    % ---------------------------
    % Método das Diferenças Concentradas
    % ---------------------------
    y_dc = zeros(1, N+1);
    y_dc(1) = y0;
    if N >= 1
        y_dc(2) = y_dc(1) + h * r * y_dc(1); % mesmo que Euler no primeiro passo
    end
    for n = 2:N
        y_dc(n+1) = 2*h*r*y_dc(n) + y_dc(n-1);
    end
    hplot = plot(t, y_dc, 's--', 'LineWidth', 1.2);
    handles(end+1) = hplot;
    labels{end+1} = sprintf('Dif. Conc. (h = %.2f)', h);
end

% ---------------------------
% Solução analítica
% ---------------------------
t_exact = linspace(t0, tf, 200);
y_exact = y0 * exp(r * t_exact);
hplot = plot(t_exact, y_exact, 'r-', 'LineWidth', 1.8);
handles(end+1) = hplot;
labels{end+1} = 'Solução Analítica';

% Ajustes do gráfico
xlabel('t');
ylabel('y(t)');
legend(handles, labels, 'Location', 'northwest');
title('Comparação: Euler vs Diferenças Concentradas');
grid on;

