from fractions import Fraction

def tree_cheeger(n, edges):
    if n < 2 or len(edges) != n - 1:
        raise ValueError("need a tree with at least 2 vertices (n - 1 edges)")

    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    parent = [-1] * n
    order = [0]
    seen = [False] * n
    seen[0] = True
    for u in order:
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                order.append(w)
    if len(order) != n:
        raise ValueError("the graph is not connected, so it is not a tree")

    size = [1] * n
    for v in reversed(order[1:]):
        size[parent[v]] += size[v]

    # most balanced edge (parent[v], v) where phi = 1 / (n - 1 - delta)
    best = min(order[1:], key=lambda v: abs(n - 2 * size[v]))
    delta = abs(n - 2 * size[best])
    phi = Fraction(1, n - 1 - delta)

    sub = set()
    stack = [best]
    while stack:
        u = stack.pop()
        sub.add(u)
        stack.extend(w for w in adj[u] if w != parent[u])
    side = sub if 2 * size[best] <= n else set(range(n)) - sub
    return phi, sorted(side)