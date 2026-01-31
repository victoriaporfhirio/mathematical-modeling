% Implementação do método da Secante no Octave
clear; clc; close all;

% Definir a função
f = @(x) 3*x - exp(x);  % Função escolhida

% Entrada dos valores iniciais
x0 = input('Digite o valor inicial x0: ');
x1 = input('Digite o valor inicial x1: ');
tol = input('Digite a tolerância desejada: ');
max_iter = input('Digite o número máximo de iterações: ');

fprintf('\nIteração\t\tx0\t\tx1\t\tx2\t\tf(x2)\n');

for i = 1:max_iter
    fx0 = f(x0);
    fx1 = f(x1);

    % Verificar divisão por zero
    if fx1 - fx0 == 0
        fprintf('Divisão por zero detectada. Encerrando.\n');
        break;
    end

    % Fórmula da secante
    x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    fx2 = f(x2);

    % Mostrar iteração
    fprintf('%d\t\t%.6f\t%.6f\t%.6f\t%.6f\n', i, x0, x1, x2, fx2);

    % Critério de parada
    if abs(fx2) < tol || abs(x2 - x1) < tol
        fprintf('\nRaiz aproximada encontrada: %.6f\n', x2);
        break;
    end

    % Atualizar valores para próxima iteração
    x0 = x1;
    x1 = x2;
end

% Comparação com fzero, se intervalo for válido
if f(x0)*f(x1) < 0
    raiz_fzero = fzero(f, [x0 x1]);
    fprintf('\nComparação com fzero: raiz = %.6f\n', raiz_fzero);
else
    fprintf('\nfzero não pôde ser aplicado: f(x0) e f(x1) têm o mesmo sinal.\n');
end

