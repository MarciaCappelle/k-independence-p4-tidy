"""
The dynamic programming algorithm over the labelled tree T(G):
Algorithm 1 (join), Algorithm 2 (main), and the closed formulas of Table 1.

Every node of T(G) is given the vector A[1..k] with A[j] = alpha_j of the
graph it represents.
"""

from p4tidy_tree import labelled_tree

# alpha_j of the exceptional graphs (Observation 20)
EXCEPTIONAL = {
    "P5":    lambda j: 3 if j == 1 else (4 if j == 2 else 5),
    "P5bar": lambda j: 2 if j == 1 else (3 if j == 2 else (4 if j == 3 else 5)),
    "C5":    lambda j: 2 if j == 1 else (3 if j == 2 else 5),
}


def join_vector(A, B, k):
    """Algorithm 1: alpha_j of a join, for j = 1..k."""
    out = [0] * (k + 1)
    for j in range(1, k + 1):
        if j == 1:
            out[1] = max(A[1], B[1])
            continue
        best = max(A[j], B[j])
        for x in range(1, j):
            for y in range(1, j):
                val = min(A[x], j - y) + min(B[y], j - x)
                if val > best:
                    best = val
        out[j] = best
    return out


def spider_formula(label, j, a):
    """Table 1.  `a[i]` is alpha_i(G[R]) (zero when the head is empty)."""
    m, r = label["m"], label["r"]
    thick, side, kind = label["thick"], label["side"], label["kind"]
    aj, a1 = a[j], a[1]

    if side is None:                                   # spider
        if not thick:
            return m + max(aj, min(j - 1, m + r))
        if m <= j and r < j - m + 1:
            return min(j + 1, 2 * m + r)
        return m + aj

    if not thick:                                      # thin quasi-spider
        if side == "X" and kind == "false":
            if j == 1:
                return m + 1 + a1
            return m + 1 + max(aj, min(j - 1, r + m - 1), min(j - 2, r + m))
        if side == "X" and kind == "true":
            if j == 1:
                return m + a1
            return m + 1 + max(aj, min(j - 1, r + m - 1), min(j - 2, r + m))
        if side == "C" and kind == "true":
            return m + max(aj, min(j - 1, r + m + 1))
        if side == "C" and kind == "false":
            if j <= 2:
                return m + max(aj, 1)
            return m + max(aj, min(j - 1, r + m + 1), min(j, r + 2))

    else:                                              # thick quasi-spider
        if side == "X" and kind == "false":
            if j < m:
                return m + 1 + aj
            return m + 1 + max(aj, min(j - m, r + m), min(j - m + 1, r + 1))
        if side == "X" and kind == "true":
            if j == 1:
                return m + a1
            if j < m:
                return m + 1 + aj
            return m + 1 + max(aj, min(j - m, r + m), min(j - m + 1, r + 1))
        if side == "C" and kind == "true":
            if j < m:
                return m + aj
            return m + max(aj, min(j - m + 1, r + m + 1))
        if side == "C" and kind == "false":
            if j <= 2:
                return max(m + aj, 3)
            if j < m:
                return max(m + aj, 3 + min(j - 1, m - 1 + a[j - 2]))
            return max(m + aj,
                       m + min(j - m + 1, r + m + 1),
                       m + 2 + min(r, j - m),
                       3 + min(j - 1, m - 1 + a[j - 2]))

    raise AssertionError("unreachable")


def dp_vector(T, k):
    """alpha_j for j = 1..k at node T, computed bottom-up."""
    if T.kind == "leaf":
        return [0] + [1] * k

    if T.kind == "union":
        out = [0] * (k + 1)
        for c in T.children:
            v = dp_vector(c, k)
            for j in range(1, k + 1):
                out[j] += v[j]
        return out

    if T.kind == "join":
        A = dp_vector(T.children[0], k)
        B = dp_vector(T.children[1], k)
        return join_vector(A, B, k)

    if T.kind == "exceptional":
        f = EXCEPTIONAL[T.label["type"]]
        return [0] + [f(j) for j in range(1, k + 1)]

    if T.kind == "spider":
        a = [0] * (k + 1)
        if T.children:
            a = dp_vector(T.children[0], k)
        return [0] + [spider_formula(T.label, j, a) for j in range(1, k + 1)]

    raise AssertionError(f"unknown node kind {T.kind}")


def alpha_k_algorithm(n, adj, k, premark=True):
    """Algorithm 2: line 1 is the shortcut k >= n(G), then the DP."""
    if k >= n:                       # Observation 4
        return n
    T = labelled_tree(n, adj, premark=premark)
    return dp_vector(T, k)[k]


def alpha_vector_algorithm(n, adj, kmax, premark=True):
    """alpha_j for j = 1..kmax, each computed through Algorithm 2."""
    return [0] + [alpha_k_algorithm(n, adj, j, premark) for j in range(1, kmax + 1)]
