from itertools import combinations
def cheeger_conductance(vertices, edges):
    vertices = list(vertices)
    n = len(vertices)

    if n < 2:
        raise ValueError("The graph needs at least 2 vertices.")

    if len(set(vertices)) != n:
        raise ValueError("Vertices must be unique.")

    vertex_set = set(vertices)
    edge_set = set()

    # Validate and normalize undirected edges.
    for u, v in edges:
        if u not in vertex_set or v not in vertex_set:
            raise ValueError("Edge endpoint is not a vertex.")
        if u == v:
            raise ValueError("Self-loops are not supported.")
        edge_set.add(frozenset((u, v)))

    # Calculate vertex degrees.
    degrees = {v: 0 for v in vertices}
    for edge in edge_set:
        u, v = tuple(edge)
        degrees[u] += 1
        degrees[v] += 1
    if any(degree == 0 for degree in degrees.values()):
        raise ValueError("Graph must have no isolated vertices.")
    total_volume = sum(degrees.values())
    half_volume = total_volume / 2

    best_ratio = float("inf")
    best_subset = None
    best_boundary = None
    best_volume = None

    # Enumerate every nonempty subset.
    for k in range(1, n + 1):
        for subset_tuple in combinations(vertices, k):
            subset = set(subset_tuple)

            # Calculate the sum of degrees in the subset.
            volume = sum(degrees[v] for v in subset)

            if volume > half_volume:
                continue

            # Count edges with exactly one endpoint in S.
            boundary = sum(
                1 for edge in edge_set
                if sum(v in subset for v in edge) == 1
            )

            ratio = boundary / volume

            if ratio < best_ratio:
                best_ratio = ratio
                best_subset = subset.copy()
                best_boundary = boundary
                best_volume = volume

    print("=== Brute-Force Conductance Analysis ===")
    print(f"Number of vertices: {n}")
    print(f"Number of edges: {len(edge_set)}")
    print(f"Total graph volume: {total_volume}")
    print(f"Minimizing subset: {best_subset}")
    print(f"Boundary edges: {best_boundary}")
    print(f"Subset volume: {best_volume}")
    print(f"Exact volume-based conductance: {best_ratio:.6f}")