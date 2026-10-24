"""
Modular decomposition and the labelled decomposition tree T(G) of Lemma 21.

The labelling procedure follows the proof of Lemma 21 literally, including the
computation of the degrees of the representative graph Q_N by scanning the
adjacency list of one representative vertex per small child and counting
distinct children with a marking array.

The flag `premark` controls the one step under discussion:

  premark=True   the child containing the representative v is marked before
                 the scan, so the edge v v' internal to a true-twin child is
                 not counted as a neighbouring child of Q_N  (correct);
  premark=False  that step is omitted, reproducing the behaviour of the
                 procedure as stated without the intra-module exclusion.
"""

from p4tidy_core import bits, popcount, components, cocomponents

# --------------------------------------------------------------------------
# modular decomposition
# --------------------------------------------------------------------------


class MDNode:
    __slots__ = ("kind", "mask", "children", "parent")

    def __init__(self, kind, mask):
        self.kind = kind          # 'leaf' | 'parallel' | 'series' | 'prime'
        self.mask = mask
        self.children = []
        self.parent = None

    def __repr__(self):
        return f"<{self.kind} {sorted(bits(self.mask))}>"


def is_module(adj, universe, S):
    """S is a module of the subgraph induced by `universe`."""
    outside = universe & ~S
    for v in bits(outside):
        inter = adj[v] & S
        if inter != 0 and inter != S:
            return False
    return True


def maximal_proper_modules(adj, universe):
    """All inclusion-maximal modules that are proper, nonempty subsets."""
    mods = []
    sub = universe
    while True:
        if sub != 0 and sub != universe and is_module(adj, universe, sub):
            mods.append(sub)
        if sub == 0:
            break
        sub = (sub - 1) & universe
    maximal = []
    for m in mods:
        if not any(o != m and (m & o) == m for o in mods):
            maximal.append(m)
    return maximal


def md_tree(adj, universe):
    """Modular decomposition tree of the subgraph induced by `universe`."""
    if popcount(universe) == 1:
        return MDNode("leaf", universe)
    comps = components(adj, universe)
    if len(comps) > 1:
        node = MDNode("parallel", universe)
        parts = comps
    else:
        cocomps = cocomponents(adj, universe)
        if len(cocomps) > 1:
            node = MDNode("series", universe)
            parts = cocomps
        else:
            node = MDNode("prime", universe)
            parts = maximal_proper_modules(adj, universe)
            covered = 0
            for p in parts:
                assert covered & p == 0, "maximal modules overlap"
                covered |= p
            assert covered == universe, "maximal modules do not cover"
    for p in parts:
        child = md_tree(adj, p)
        child.parent = node
        node.children.append(child)
    return node


# --------------------------------------------------------------------------
# the map s(.) of Lemma 21
# --------------------------------------------------------------------------


def build_s_map(root):
    """s(y) = the small child of a prime node containing y, when it exists."""
    s = {}

    def walk(node):
        par = node.parent
        if node.kind == "leaf" and par is not None and par.kind == "prime":
            s[next(bits(node.mask))] = node
        elif (popcount(node.mask) == 2 and par is not None
              and par.kind == "prime"):
            for v in bits(node.mask):
                s[v] = node
        for ch in node.children:
            walk(ch)

    walk(root)
    return s


def representative_degrees(node, adj, s_map, premark=True):
    """Degrees of the representative graph Q_N of the prime node `node`."""
    small = [c for c in node.children if popcount(c.mask) <= 2]
    big = [c for c in node.children if popcount(c.mask) > 2]
    assert len(big) <= 1, "more than one non-small child of a prime node"
    HN = big[0] if big else None

    deg = {}
    hn_deg = 0
    for c in small:
        v = next(bits(c.mask))                      # representative vertex
        marked = set()
        if premark:
            marked.add(id(s_map[v]))                # exclude v's own child
        count = 0
        for y in bits(adj[v]):
            sy = s_map.get(y)
            if sy is not None and sy.parent is node:
                target = sy
            elif (node.mask >> y) & 1:
                target = HN
            else:
                continue                            # outside M(N)
            if target is None:
                continue
            if id(target) not in marked:
                marked.add(id(target))
                count += 1
                if target is HN:
                    hn_deg += 1
        deg[id(c)] = count
    if HN is not None:
        deg[id(HN)] = hn_deg
    return deg


# --------------------------------------------------------------------------
# the labelled tree T(G)
# --------------------------------------------------------------------------


class TNode:
    __slots__ = ("kind", "n", "children", "label")

    def __init__(self, kind, n, label=None):
        self.kind = kind      # 'leaf'|'union'|'join'|'exceptional'|'spider'
        self.n = n
        self.children = []
        self.label = label or {}

    def __repr__(self):
        return f"<T {self.kind} n={self.n} {self.label}>"


DEGSEQ = {
    (2, 2, 2, 2, 2): "C5",
    (1, 1, 2, 2, 2): "P5",
    (2, 2, 2, 3, 3): "P5bar",
    (1, 1, 2, 3, 3): "bull",
}


class LabellingError(Exception):
    pass


def build_T(md, adj, s_map, premark=True):
    """Turn the modular decomposition tree into the labelled tree T(G)."""
    if md.kind == "leaf":
        return TNode("leaf", 1)

    if md.kind == "parallel":
        node = TNode("union", popcount(md.mask))
        node.children = [build_T(c, adj, s_map, premark) for c in md.children]
        return node

    if md.kind == "series":
        # replaced by a left-deep binary tree of join nodes
        kids = [build_T(c, adj, s_map, premark) for c in md.children]
        cur = kids[0]
        for nxt in kids[1:]:
            j = TNode("join", cur.n + nxt.n)
            j.children = [cur, nxt]
            cur = j
        return cur

    # prime node
    deg = representative_degrees(md, adj, s_map, premark)
    children = md.children
    p = len(children)

    if p == 5 and all(c.kind == "leaf" for c in children):
        seq = tuple(sorted(deg[id(c)] for c in children))
        name = DEGSEQ.get(seq)
        if name is None:
            raise LabellingError(f"unknown 5-vertex quotient, degrees {seq}")
        if name != "bull":
            return TNode("exceptional", 5, {"type": name})
        # the bull is the five-vertex thin spider; fall through

    t = p // 2
    h = p - 2 * t
    if h not in (0, 1) or t < 2:
        raise LabellingError(f"quotient is not a prime spider: p={p}")

    degs = {id(c): deg[id(c)] for c in children}
    thin = any(d == 1 for d in degs.values())

    if thin:
        dX, dHead, dC = 1, t, t + h
    else:
        dX, dHead, dC = t - 1, t, 2 * t - 2 + h

    Xs, Cs, heads = [], [], []
    for c in children:
        d = degs[id(c)]
        if d == dX:
            Xs.append(c)
        elif h == 1 and d == dHead and len(heads) == 0:
            heads.append(c)
        elif d == dC:
            Cs.append(c)
        else:
            raise LabellingError(
                f"degree {d} fits no position (thin={thin}, t={t}, h={h})")

    if len(Xs) != t or len(Cs) != t or len(heads) != h:
        raise LabellingError(
            f"bad partition sizes |X|={len(Xs)} |C|={len(Cs)} "
            f"|head|={len(heads)} (t={t}, h={h})")

    side = kind = None
    for c in Xs:
        if popcount(c.mask) == 2:
            side, kind = "X", ("true" if c.kind == "series" else "false")
    for c in Cs:
        if popcount(c.mask) == 2:
            if side is not None:
                raise LabellingError("two twin pairs in one quotient")
            side, kind = "C", ("true" if c.kind == "series" else "false")

    head = heads[0] if heads else None
    r = popcount(head.mask) if head is not None else 0
    node = TNode("spider", popcount(md.mask), {
        "thick": not thin, "m": t, "r": r, "side": side, "kind": kind,
    })
    if head is not None:
        node.children = [build_T(head, adj, s_map, premark)]
    return node


def labelled_tree(n, adj, premark=True):
    universe = (1 << n) - 1
    md = md_tree(adj, universe)
    s_map = build_s_map(md)
    return build_T(md, adj, s_map, premark)


def count_nodes(T):
    return 1 + sum(count_nodes(c) for c in T.children)
