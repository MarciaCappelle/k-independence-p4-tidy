"""
Core graph utilities, brute-force k-independence, a direct P_4-tidy test,
and constructors/generators for P_4-tidy graphs.

Graphs are represented as (n, adj) where adj is a list of n bitmasks:
bit j of adj[i] is set iff ij is an edge.  Vertex sets are bitmasks too.

This module is self-contained and has no third-party dependencies.
"""

from itertools import combinations

# --------------------------------------------------------------------------
# basic utilities
# --------------------------------------------------------------------------


def bits(mask):
    """Iterate over the indices set in mask."""
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def popcount(mask):
    return bin(mask).count("1")


def induced(n, adj, mask):
    """Induced subgraph on `mask`, renumbered to 0..|mask|-1.

    Returns (m, adj2, verts) where verts[i] is the original vertex of i.
    """
    verts = list(bits(mask))
    index = {v: i for i, v in enumerate(verts)}
    adj2 = []
    for v in verts:
        row = 0
        for u in bits(adj[v] & mask):
            row |= 1 << index[u]
        adj2.append(row)
    return len(verts), adj2, verts


def complement(n, adj):
    full = (1 << n) - 1
    return n, [(~adj[v]) & full & ~(1 << v) for v in range(n)]


def disjoint_union(g1, g2):
    n1, a1 = g1
    n2, a2 = g2
    return n1 + n2, list(a1) + [row << n1 for row in a2]


def join(g1, g2):
    n1, a1 = g1
    n2, a2 = g2
    mask1 = (1 << n1) - 1
    mask2 = ((1 << n2) - 1) << n1
    adj = [a1[v] | mask2 for v in range(n1)]
    adj += [(a2[v] << n1) | mask1 for v in range(n2)]
    return n1 + n2, adj


def add_edge(adj, u, v):
    adj[u] |= 1 << v
    adj[v] |= 1 << u


# --------------------------------------------------------------------------
# connectivity
# --------------------------------------------------------------------------


def components(adj, mask):
    """Connected components of the subgraph induced by `mask`, as bitmasks."""
    comps = []
    todo = mask
    while todo:
        start = todo & -todo
        comp = start
        frontier = start
        while frontier:
            nxt = 0
            for v in bits(frontier):
                nxt |= adj[v] & mask & ~comp
            comp |= nxt
            frontier = nxt
        comps.append(comp)
        todo &= ~comp
    return comps


def cocomponents(adj, mask):
    """Connected components of the complement of the subgraph induced by mask."""
    comps = []
    todo = mask
    while todo:
        start = todo & -todo
        comp = start
        frontier = start
        while frontier:
            nxt = 0
            for v in bits(frontier):
                nxt |= (~adj[v]) & mask & ~comp & ~(1 << v)
            comp |= nxt
            frontier = nxt
        comps.append(comp)
        todo &= ~comp
    return comps


# --------------------------------------------------------------------------
# brute-force k-independence
# --------------------------------------------------------------------------


def alpha_all(n, adj, kmax=None):
    """Brute force.  Returns a list A with A[k] = alpha_k(G) for k = 1..kmax.

    A single pass over all 2^n subsets computes every k at once: for a subset
    S let d(S) be the maximum degree of G[S]; then S is k-independent iff
    d(S) <= k-1.
    """
    if kmax is None:
        kmax = n
    best = [0] * (n + 2)  # best[d] = max |S| with max-degree exactly <= d
    for s in range(1 << n):
        d = 0
        for v in bits(s):
            dv = popcount(adj[v] & s)
            if dv > d:
                d = dv
        size = popcount(s)
        if size > best[d]:
            best[d] = size
    # prefix maxima: alpha_k = max over d <= k-1
    run = 0
    cum = []
    for d in range(n + 2):
        run = max(run, best[d])
        cum.append(run)
    return [0] + [cum[min(k - 1, n + 1)] for k in range(1, kmax + 1)]


def alpha_k_bruteforce(n, adj, k):
    return alpha_all(n, adj, k)[k]


# --------------------------------------------------------------------------
# direct P_4-tidy test (from the definition)
# --------------------------------------------------------------------------


def induces_p4(adj, quad):
    """True iff the four given vertices induce a P_4."""
    degs = []
    edges = 0
    for v in quad:
        d = 0
        for u in quad:
            if u != v and (adj[v] >> u) & 1:
                d += 1
        degs.append(d)
        edges += d
    return edges // 2 == 3 and sorted(degs) == [1, 1, 2, 2]


def is_p4_tidy(n, adj):
    """Direct check: every vertex set A inducing a P_4 has at most one partner,
    i.e. at most one v not in A such that A + v induces at least two P_4's."""
    quads = [q for q in combinations(range(n), 4) if induces_p4(adj, q)]
    for q in quads:
        qs = set(q)
        partners = 0
        for v in range(n):
            if v in qs:
                continue
            five = list(q) + [v]
            count = sum(1 for sub in combinations(five, 4) if induces_p4(adj, sub))
            if count >= 2:
                partners += 1
                if partners > 1:
                    return False
    return True


# --------------------------------------------------------------------------
# constructors
# --------------------------------------------------------------------------


def K1():
    return 1, [0]


def path(n):
    adj = [0] * n
    for i in range(n - 1):
        add_edge(adj, i, i + 1)
    return n, adj


def cycle(n):
    nn, adj = path(n)
    add_edge(adj, 0, n - 1)
    return nn, adj


def P5():
    return path(5)


def P5bar():
    return complement(*path(5))


def C5():
    return cycle(5)


def spider(t, thick=False, head=None, twin_side=None, twin_kind=None):
    """Build a spider or quasi-spider of weight t with the given head.

    Layout: X = 0..t-1, C = t..2t-1, then the twin copy (if any), then head.
    twin_side in {None,'X','C'}; twin_kind in {'true','false'}.
    The duplicated vertex is always the one with index s = t-1.
    """
    assert t >= 2
    size = 2 * t
    adj = [0] * size
    X = list(range(t))
    C = list(range(t, 2 * t))
    for i in range(t):
        for j in range(i + 1, t):
            add_edge(adj, C[i], C[j])
    for i in range(t):
        if thick:
            for j in range(t):
                if j != i:
                    add_edge(adj, X[i], C[j])
        else:
            add_edge(adj, X[i], C[i])

    g = (size, adj)

    # duplicate vertex s into a twin pair
    if twin_side is not None:
        s = (X if twin_side == "X" else C)[t - 1]
        n0, a0 = g
        new = n0
        a0 = list(a0) + [0]
        n0 += 1
        for u in bits(a0[s]):
            add_edge(a0, new, u)
        if twin_kind == "true":
            add_edge(a0, new, s)
        g = (n0, a0)
        C = list(C)
        X = list(X)

    # attach the head: complete to C, anticomplete to X (and to the twin copy
    # when the twin pair sits in X)
    if head is not None:
        nh, ah = head
        n0, a0 = g
        a0 = list(a0) + [0] * nh
        for v in range(nh):
            for u in bits(ah[v]):
                add_edge(a0, n0 + v, n0 + u)
        cmask = 0
        for c in C:
            cmask |= 1 << c
        if twin_side == "C":
            cmask |= 1 << (2 * t)  # the twin copy also belongs to C
        for v in range(nh):
            for c in bits(cmask):
                add_edge(a0, n0 + v, c)
        g = (n0 + nh, a0)

    return g


def canonical_key(n, adj):
    return (n, tuple(adj))


def all_graphs_upto(nmax):
    """One representative of every isomorphism class of simple graphs with
    1..nmax vertices (brute force; nmax is meant to be small)."""
    from itertools import permutations
    out = []
    for n in range(1, nmax + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        seen = set()
        for mask in range(1 << len(pairs)):
            adj = [0] * n
            for b, (i, j) in enumerate(pairs):
                if (mask >> b) & 1:
                    add_edge(adj, i, j)
            canon = min(
                tuple(sorted(
                    (min(p[i], p[j]), max(p[i], p[j]))
                    for i in range(n) for j in bits(adj[i]) if i < j))
                for p in permutations(range(n)))
            if canon in seen:
                continue
            seen.add(canon)
            out.append((n, adj))
    return out


def generate_corpus(max_n=13, verbose=False):
    """Generate a corpus of P_4-tidy graphs built from the structure theorem,
    covering nonempty heads, true/false twin modules on both sides, and nested
    union/join operations."""
    base = [K1(), P5(), P5bar(), C5(), path(2), path(3), path(4)]
    corpus = []
    seen = set()

    def push(g, tag):
        n, adj = g
        if n > max_n:
            return
        key = canonical_key(n, adj)
        if key in seen:
            return
        seen.add(key)
        corpus.append((g, tag))

    for g in base:
        push(g, "base")

    # spiders and quasi-spiders, with and without heads, nested heads included
    # Heads are of two kinds: every simple graph on at most four vertices,
    # up to isomorphism, and structured heads that are themselves spiders,
    # quasi-spiders, unions and joins, so that decompositions are nested
    # several levels deep.
    heads = [None]
    heads += all_graphs_upto(4)
    heads += [P5(), C5(), P5bar(), disjoint_union(P5(), K1()),
              spider(2, False),                       # a P_4 as head
              spider(2, False, None, "X", "true"),    # quasi-spider as head
              spider(2, False, None, "C", "false"),
              spider(3, True),                        # thick spider as head
              join(path(2), complement(*path(2)))]
    for t in (2, 3, 4):
        for thick in (False, True):
            for head in heads:
                push(spider(t, thick, head), f"spider t={t} thick={thick}")
                for side in ("X", "C"):
                    for kind in ("true", "false"):
                        push(spider(t, thick, head, side, kind),
                             f"quasi t={t} thick={thick} {side}^{kind}")

    # nested union / join over everything built so far
    seeds = [g for g, _ in corpus if g[0] <= 7]
    for g1 in seeds:
        for g2 in seeds:
            push(disjoint_union(g1, g2), "union")
            push(join(g1, g2), "join")

    if verbose:
        print(f"corpus: {len(corpus)} graphs")
    return corpus


if __name__ == "__main__":
    corpus = generate_corpus(max_n=12, verbose=True)
    bad = [(g, tag) for g, tag in corpus if not is_p4_tidy(*g)]
    print("non-P4-tidy in corpus:", len(bad))
    for g, tag in bad[:5]:
        print("  ", tag, g[0])
    sizes = {}
    for g, tag in corpus:
        sizes[g[0]] = sizes.get(g[0], 0) + 1
    print("sizes:", dict(sorted(sizes.items())))
