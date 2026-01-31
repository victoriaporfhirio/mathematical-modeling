# importações de bibliotecas
import itertools         # utilitários para combinações
import numpy as np       # numpy para manipulação de matrizes/arrays
import pulp              # PuLP: modelagem de PL/PI e interface com CBC
import math              # utilitários matemáticos
from copy import deepcopy # copia profunda de listas/dicionários (ramificar)

# Função para o usuário DIGITAR a matriz
def ler_matriz(): # pede ao usuário o número de cidades (dimensão da matriz)
    print("Digite o valor de n (número de cidades): ")
    n = int(input())                # lê "n" como inteiro

 # instruções de como digitar a matriz (exemplo)
    print("\nDigite a matriz de distâncias linha por linha:")
    #print("Ex: para n=3 digite:")
    #print("0 4 9")
   # print("4 0 2")
   # print("9 2 0\n")

    matriz = []                     # lista temporária para armazenar linhas
    for i in range(n):              # para cada linha i = 0,..,n-1
        linha = input().strip().split()   # lê a linha, remove espaços das pontas e separa por espaços
        linha = list(map(float, linha))   # converte cada token para float (permite custos decimais)
        matriz.append(linha)              # adiciona a linha lida à lista

    return np.array(matriz)         # converte para numpy.array (matriz) e retorna

# Função que resolve o LP relaxado com cortes SEC
def solve_lp_relax(dist, cuts=None, fixed_vars=None, msg=False):
# dist: matriz de custos (nxn)
# cuts: lista de conjuntos S (cada S é um set de índices) para restrições SEC
# fixed_vars: dicionário {(i,j): 0/1} com variáveis fixadas (ramificações)
# msg: se True, mostra mensagens do solver

    if cuts is None:
        cuts = []                    # evita mutação do argumento padrão
    if fixed_vars is None:
        fixed_vars = {}              # idem

    n = len(dist)                   # dimensão da matriz (n cidades)
    prob = pulp.LpProblem("TSP_relax", pulp.LpMinimize)  # cria problema LP com objetivo minimizar

    # cria variáveis x[i][j] contínuas no intervalo [0,1] (relaxação)
    x = pulp.LpVariable.dicts("x", (range(n), range(n)), lowBound=0, upBound=1, cat=pulp.LpContinuous)

    # função objetivo: minimizar soma dos custos * x_ij
    prob += pulp.lpSum(dist[i][j] * x[i][j] for i in range(n) for j in range(n))

    # proíbe self-loops: x[i][i] == 0 para todo i
    for i in range(n):
        prob += x[i][i] == 0

    # restrição: saída de cada cidade = 1 (exatamente uma aresta saindo)
    for i in range(n):
        prob += pulp.lpSum(x[i][j] for j in range(n) if j != i) == 1

    # restrição: entrada em cada cidade = 1 (exatamente uma aresta entrando)
    for j in range(n):
        prob += pulp.lpSum(x[i][j] for i in range(n) if i != j) == 1

    # adicionar cortes SEC (Subtour Elimination Constraints) passados em cuts
    for S in cuts:
        prob += pulp.lpSum(x[i][j] for i in S for j in S if i != j) <= len(S) - 1

    # aplicar fixings vindos do branching (forçar x_ij = 0 ou 1)
    for (i, j), val in fixed_vars.items():
        if val == 1:
            prob += x[i][j] == 1
        else:
            prob += x[i][j] == 0

    # resolve o LP usando CBC via PuLP; msg controla se imprime logs
    prob.solve(pulp.PULP_CBC_CMD(msg=1 if msg else 0))

    # captura o status (ex.: "Optimal") e valor objetivo (se ótimo)
    status = pulp.LpStatus[prob.status]
    obj = pulp.value(prob.objective) if prob.status == 1 else float('inf')

    # extrai a matriz solução com valores das variáveis x[i][j]
    sol = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            v = x[i][j].value()           # valor numérico da variável (pode ser None)
            sol[i][j] = 0.0 if v is None else float(v)

    return sol, obj, status            # retorna matriz de valores, custo e status

# Detectar subtours
def find_subtours_from_solution(sol): # sol: matriz (n x n) de valores x_ij (contínuos)
    n = sol.shape[0]
    unvisited = set(range(n))          # nós/índices ainda não explorados
    subtours = []                      # lista onde guardaremos os subtours encontrados

    # enquanto houver nós não visitados, constrói-se um caminho/ciclo
    while unvisited:
        start = unvisited.pop()        # pega um nó inicial
        tour = [start]                 # caminho temporário
        cur = start                    # nó corrente
        while True:
            # escolhe j que maximiza sol[cur][j], ignorando self-loop (usando -1 para evitar escolher cur)
            j = int(np.argmax([sol[cur][k] if k != cur else -1 for k in range(n)]))
            if j in tour:
                # fechou um ciclo: encontra índice de primeira ocorrência de j no caminho
                idx = tour.index(j)
                cycle = tour[idx:]         # extrai o ciclo (subtour)
                subtours.append(cycle)     # adiciona o subtour à lista
                # remove os vértices do ciclo do conjunto de não visitados
                for v in cycle:
                    if v in unvisited:
                        unvisited.remove(v)
                break
            else:
                # continua seguindo o arco mais forte
                tour.append(j)
                if j in unvisited:
                    unvisited.remove(j)
                cur = j

    return subtours                     # retorna lista de subtours (cada um é lista de índices)

# Checar integridade
def is_integer_solution(sol, tol=1e-6):
    # verifica se todos os valores da matriz sol são essencialmente 0 ou 1
    n = sol.shape[0]
    for i in range(n):
        for j in range(n):
            if abs(sol[i][j] - round(sol[i][j])) > tol:
                return False
    return True

# Extrair tour de solução inteira
def extract_tour_from_integral_sol(sol):
    # constrói mapa de sucessores succ[i] = j onde x[i][j] == 1
    n = sol.shape[0]
    succ = {}
    for i in range(n):
        for j in range(n):
            if round(sol[i][j]) == 1:   # se x[i][j] é 1 (inteiro)
                succ[i] = j
                break                    # cada linha tem exatamente uma saída

    # reconstruir a sequência do tour começando do vértice 0
    tour = []
    cur = 0
    for _ in range(n):
        tour.append(cur)
        cur = succ[cur]
    return tour

def tour_cost(tour, dist):
    # calcula custo total do tour (soma dist[i][j] ao percorrer o ciclo)
    cost = 0
    for k in range(len(tour)):
        i = tour[k]
        j = tour[(k + 1) % len(tour)]  # % fecha o ciclo de volta ao início
        cost += dist[i][j]
    return cost

# BRANCH-AND-CUT COM SEC
def branch_and_cut_sec(dist, verbose=False, max_nodes=5000):
    # dist: matriz de custos (numpy array)
    # verbose: se True imprime logs de progresso
    # max_nodes: limite de nós a explorar (proteção contra explosão)
    n = len(dist)
    initial_state = ([], {})            # estado inicial: (lista de cortes, dicionário de fixings)
    stack = [initial_state]             # pilha para DFS (branch-and-bound)
    incumbent_solution = None           # melhor solução inteira encontrada
    incumbent_cost = float('inf')       # custo do incumbent
    nodes_explored = 0                  # contador de nós processados

    # enquanto houver estados na pilha
    while stack:
        nodes_explored += 1
        cuts, fixed = stack.pop()      # retira o último estado empilhado (DFS)

        # resolve a relaxação LP no estado atual (com cortes e fixings)
        sol, obj, status = solve_lp_relax(dist, cuts=cuts, fixed_vars=fixed, msg=False)
        if status != "Optimal":
            # se LP não foi resolvido ótimo (ex.: infeasible), pula o ramo
            continue

        # poda por bound inferior: se obj LP >= melhor inteiro conhecido, não prossegue
        if obj >= incumbent_cost:
            continue

        # detecta subtours na solução LP
        subtours = find_subtours_from_solution(sol)

        if len(subtours) > 1:
            # se há mais de um subtour, cria cortes SEC correspondentes e empilha o estado com esses cortes
            new_cuts = deepcopy(cuts)   # copia a lista de cortes atual
            added = 0
            for S in subtours:
                Sset = set(S)
                # adiciona apenas cortes novos e válidos (1 <= |S| < n)
                if Sset not in new_cuts and 1 <= len(Sset) < n:
                    new_cuts.append(Sset)
                    added += 1
            if added > 0:
                # empilha o mesmo estado com cortes adicionais (faça cutting antes do branching)
                stack.append((new_cuts, deepcopy(fixed)))
                continue    # processa próximo nó da pilha

        # se a solução LP é inteira (todos 0/1), atualiza incumbent se melhor
        if is_integer_solution(sol):
            tour = extract_tour_from_integral_sol(sol)
            cost = tour_cost(tour, dist)
            if cost < incumbent_cost:
                incumbent_cost = cost
                incumbent_solution = (tour, sol.copy())
            continue

        # solução fracionária sem cortes novos: escolher variável fracionária para branch
        frac_vars = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    v = sol[i][j]
                    # considera fracionárias aquelas estritamente entre 0 e 1 (com tolerância)
                    if 1e-6 < v < 1 - 1e-6:
                        frac_vars.append(((i, j), abs(v - 0.5)))  # guarda também distância de 0.5

        if not frac_vars:
            # nenhuma variável fracionária encontrada (caso raro), pula
            continue

        # escolhe variável mais fracionária (mais próxima de 0.5) como heurística
        (i_sel, j_sel), _ = min(frac_vars, key=lambda x: x[1])

        # cria dois ramos: fixar x_ij = 1 e fixar x_ij = 0
        fixed1 = deepcopy(fixed)
        fixed1[(i_sel, j_sel)] = 1
        fixed2 = deepcopy(fixed)
        fixed2[(i_sel, j_sel)] = 0

        # empilha primeiro o ramo x=0 e depois o ramo x=1 para que x=1 seja processado antes (LIFO)
        stack.append((deepcopy(cuts), fixed2))
        stack.append((deepcopy(cuts), fixed1))

    # fim do while: retorna a melhor solução inteira encontrada (se existir)
    if incumbent_solution is None:
        return None, float('inf'), None

    tour, solmat = incumbent_solution
    return tour, incumbent_cost, solmat

# EXECUÇÃO INTERATIVA
print("=== Modo interativo TSP com Branch-and-Cut ===")
dist = ler_matriz()
print("\nMatriz carregada:\n", dist)

# executa o Branch-and-Cut com verbose=True para acompanhar logs
tour, cost, sol = branch_and_cut_sec(dist, verbose=True)

# imprime resultado final
print("\n========================")
if tour is None:
    print("Nenhuma solução inteira encontrada.")
else:
    print("Tour ótimo:", tour)
    print("Custo:", cost)
    print("Matriz solução inteira:\n", sol.round().astype(int))
