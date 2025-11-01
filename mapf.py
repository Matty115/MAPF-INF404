from pysat.formula import WCNF
from pysat.solvers import Glucose3
import heapq
from pysat.card import CardEnc, EncType

def create_mapf_cnf(agents, grid_size, obstacles, max_time):
    wcnf = WCNF()

    def d(u, v, cost_fn=None):
        if u == v:
            return 0.0

        # movimientos (arriba, abajo, izquierda, derecha)
        moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]

        # prioridad (distancia acumulada, nodo)
        pq = [(0.0, u)]
        visited = set()
        dist = {u: 0.0}

        while pq:
            d_u, (x, y) = heapq.heappop(pq)
            if (x, y) in visited:
                continue
            visited.add((x, y))

            if (x, y) == v:
                return d_u

            for dx, dy in moves:
                nx, ny = x + dx, y + dy

                # límites y obstáculos
                if not (0 <= nx < grid_size[0] and 0 <= ny < grid_size[1]):
                    continue
                if (nx, ny) in obstacles:
                    continue

                # costo de transición
                cost = cost_fn((x, y), (nx, ny)) if cost_fn else 1.0
                nd = d_u + cost

                if nd < dist.get((nx, ny), float('inf')):
                    dist[(nx, ny)] = nd
                    heapq.heappush(pq, (nd, (nx, ny)))

        return float('inf')
    
    def feasible(a, t):
        nodes = []
        for x in range(grid_size[0]):
            for y in range(grid_size[1]):
                if d(agents[a][0], (x, y)) <= t and d((x, y), agents[a][1]) <= (max_time - t):
                    nodes.append((x, y))
        return nodes
    
    def neighbors(u):
        x, y = u
        moves = [(0,0),(1,0),(-1,0),(0,1),(0,-1)]
        valid = []
        for dx, dy in moves:
            nx, ny = x+dx, y+dy
            if 0 <= nx < grid_size[0] and 0 <= ny < grid_size[1] and (nx, ny) not in obstacles:
                valid.append((nx, ny))
        return valid

    

    # Variable encoding
    def at(a, u, t):
        return a * (grid_size[0] * grid_size[1] * (max_time + 1)) + u[0] * (
            grid_size[1] * (max_time + 1)) + u[1] * (max_time + 1) + t + 1
    
    def shift(u, v, t):
        return at(len(agents), 0, 0) + u[0] * (
            grid_size[0] * grid_size[1]**2 * (max_time + 1)) + u[1] * (
            grid_size[0] * grid_size[1] * (max_time + 1)) + v[0] * (
            grid_size[1] * (max_time + 1)) + v[1] * (max_time + 1) + t + 1
    
    def finalState(a, t):
        return shift(grid_size[0], 0, 0) + a * (max_time + 1) + t + 1
    

    # feasible nodes precomputation
    feasible_nodes = [[feasible(a, t) for t in range(max_time + 1)] for a in range(len(agents))]


    ## Restricciones duras y de cardinalidad
    for t in range(max_time):
        for x in range(grid_size[0]):
            for y in range(grid_size[1]):
                u = (x, y)
                lits = []
                for v in neighbors(u):
                    lits.append(shift(u, v, t))
                    # H9 constrain
                    wcnf.append([-shift(u, v, t), -shift(v, u, t)])
                # C1 constrain
                res = CardEnc.equals(lits=lits, bound=1, encoding=EncType.pairwise)
                wcnf.extend(res.clauses)

    for t in range(max_time + 1):
        for x in range(grid_size[0]):
            for y in range(grid_size[1]):
                u = (x, y)
                lits = []
                for a in range(len(agents)):
                    if u in feasible_nodes[a][t]:
                        lits.append(at(a, u, t))
                # C2 constrain
                if lits:
                    res = CardEnc.atmost(lits=lits, bound=1, encoding=EncType.ladder)
                    wcnf.extend(res.clauses)
                for v in neighbors(u):
                    # H10 constrain
                    wcnf.append([-shift(u, v, t), shift(v, v, t)])

    for a in range(len(agents)):
        # H5 constrain
        wcnf.append([finalState(a, max_time)])
        # H7 constrain
        wcnf.append([at(a, agents[a][0], 0)])
        # H8 constrain
        wcnf.append([at(a, agents[a][1], max_time)])
        for t in range(max_time):
            for u in feasible_nodes[a][t]:
                c = [-at(a, u, t)]
                c += [at(a, v, t + 1) for v in feasible_nodes[a][t] if v in neighbors(u)]
                # H11 constrain
                wcnf.append(c)
                for v in feasible_nodes[a][t + 1]:
                    if v in neighbors(u):
                        # H1, H2 constrains
                        wcnf.append([-at(a, u, t), -shift(u, v, t), at(a, v, t + 1)])
                        wcnf.append([-at(a, u, t), -at(a, v, t + 1), shift(u, v, t)])
                feasible_next = set(feasible_nodes[a][t + 1])
                for v in neighbors(u):
                    if v not in feasible_next:
                        # H4 constrain
                        wcnf.append([-at(a, u, t), -shift(u, v, t)])
            for v in feasible_nodes[a][t + 1]:
                c = [-at(a, v, t + 1)]
                if feasible_nodes[a][t]:
                    c += [at(a, u, t) for u in feasible_nodes[a][t] if v in neighbors(u)]
                # H3 constrain
                wcnf.append(c)
        for t in range(d(agents[a][0], agents[a][1]), max_time):
            # H6 constrains
            wcnf.append([-at(a, agents[a][1], t), -finalState(a, t + 1), finalState(a, t)])
            wcnf.append([-finalState(a, t), at(a, agents[a][1], t)])
            wcnf.append([-finalState(a, t), finalState(a, t + 1)])
        for t in range(max_time + 1):
            lit = [at(a, u, t) for u in feasible_nodes[a][t]]
            # C3 constrain
            res = CardEnc.equals(lits=lit, bound=1, encoding=EncType.pairwise)
            wcnf.extend(res.clauses)

        ## Restricciones blandas (maximizar costo)
        for a in range(len(agents)):
            for t in range(d(agents[a][0], agents[a][1]) ,max_time + 1):
                # S1 constrain
                wcnf.append([finalState(a, t)], weight=1)


    return wcnf