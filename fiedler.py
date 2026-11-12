import numpy as np


def fiedler_sweep_cheeger(vertices, edges):
    """
    Approximate graph conductance using the Fiedler vector
    and sweep cuts. Prints results without returning values.
    """
    vertices = list(vertices)
    n = len(vertices)

    if n < 2:
        raise ValueError("The graph needs at least 2 vertices.")
    if len(set(vertices)) != n:
        raise ValueError("Vertices must be unique.")

    index = {v: i for i, v in enumerate(vertices)}
    edge_set = set()

    for u, v in edges:
        if u not in index or v not in index:
            raise ValueError("Edge endpoint is not a vertex.")
        i, j = index[u], index[v]

        if i == j:
            raise ValueError("Self-loops are not supported.")

        edge_set.add((min(i, j), max(i, j)))

    # Construct adjacency matrix.
    A = np.zeros((n, n), dtype=float)

    for i, j in edge_set:
        A[i, j] = A[j, i] = 1.0

    degrees = A.sum(axis=1)

    if np.any(degrees == 0):
        raise ValueError("Graph must have no isolated vertices.")

    # Normalized graph Laplacian.
    inv_sqrt_deg = 1.0 / np.sqrt(degrees)
    L_norm = (
        np.eye(n)
        - inv_sqrt_deg[:, None]
        * A
        * inv_sqrt_deg[None, :]
    )

    eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
    lambda2 = max(0.0, float(eigenvalues[1]))

    # Fiedler sweep scores.
    fiedler = eigenvectors[:, 1]
    order = np.argsort(fiedler)

    total_volume = float(degrees.sum())
    best_ratio = float("inf")
    best_subset = None
    best_boundary = None
    best_volume = None

    # Adjacency lists for efficient boundary updates.
    adjacency = [[] for _ in range(n)]
    for i, j in edge_set:
        adjacency[i].append(j)
        adjacency[j].append(i)

    selected = set()
    volume = 0.0
    boundary = 0

    # Sweep cuts, maintaining volume <= half total volume.
    for idx in order:
        inside_neighbors = sum(
            neighbor in selected
            for neighbor in adjacency[idx]
        )

        volume = volume + degrees[idx]
        boundary = (
            boundary + degrees[idx] - 2 * inside_neighbors
        )
        selected.add(idx)

        if volume <= total_volume / 2:
            candidate = set(selected)
            candidate_volume = volume
            candidate_boundary = boundary
        else:
            candidate = set(range(n)) - selected
            candidate_volume = total_volume - volume
            candidate_boundary = boundary

        if candidate_volume <= 0:
            continue
        ratio = candidate_boundary / candidate_volume

        if ratio < best_ratio:
            best_ratio = ratio
            best_subset = {vertices[i] for i in candidate}
            best_boundary = int(candidate_boundary)
            best_volume = candidate_volume

    # Cheeger inequality bounds for conductance.
    lower_bound = lambda2 / 2
    upper_bound = np.sqrt(2 * lambda2)

    # Print the main results.
    print("=== Fiedler Sweep Cheeger Analysis ===")
    print(f"Number of vertices: {n}")
    print(f"Number of edges: {len(edge_set)}")
    print(f"Second eigenvalue (lambda_2): {lambda2:.6f}")
    print(f"Best subset: {best_subset}")
    print(f"Boundary edges: {best_boundary}")
    print(f"Subset volume: {best_volume:g}")
    print(f"Sweep conductance estimate: {best_ratio:.6f}")
    print(f"Cheeger lower bound: {lower_bound:.6f}")
    print(f"Cheeger upper bound: {upper_bound:.6f}")