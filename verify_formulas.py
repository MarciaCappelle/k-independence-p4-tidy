#!/usr/bin/env python3
"""Exhaustive consistency check for the formulas in Table 1.

The script generates all labelled simple graphs on the head R up to a
specified order, builds the two spider types and the eight quasi-spider
types, computes alpha_k exactly by exhaustive enumeration of vertex
subsets, and compares those values with the closed formulas in the
manuscript "k-independence in P4-tidy graphs".

Default parameters reproduce the verification reported in the revised
manuscript:
    m in {2,3,4}, 0 <= |R| <= 4, and 1 <= k <= |V(G)|.

No third-party packages are required.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple


Adj = Tuple[int, ...]  # adjacency bit masks


@dataclass(frozen=True)
class Structure:
    family: str              # "spider" or "quasi"
    kind: str                # "thin" or "thick"
    duplicated_side: Optional[str] = None  # None, "C", or "X"
    twin_type: Optional[str] = None        # None, "true", or "false"

    @property
    def label(self) -> str:
        if self.family == "spider":
            return f"{self.kind} spider"
        suffix = "t" if self.twin_type == "true" else "f"
        if self.duplicated_side == "C":
            return f"{self.kind} quasi (R,C^{suffix},X)"
        return f"{self.kind} quasi (R,C,X^{suffix})"


def add_edge(adj: List[int], u: int, v: int) -> None:
    adj[u] |= 1 << v
    adj[v] |= 1 << u


def all_graphs(r: int) -> Iterator[Adj]:
    """Generate every labelled simple graph on vertices 0,...,r-1."""
    pairs = list(itertools.combinations(range(r), 2))
    for mask in range(1 << len(pairs)):
        adj = [0] * r
        for j, (u, v) in enumerate(pairs):
            if mask & (1 << j):
                add_edge(adj, u, v)
        yield tuple(adj)


def spider(head: Adj, m: int, kind: str) -> Tuple[Adj, Tuple[int, ...], Tuple[int, ...]]:
    """Build a labelled thin or thick spider with head ``head``."""
    r = len(head)
    n = r + 2 * m
    adj = [0] * n

    # Copy G[R].
    for u in range(r):
        adj[u] = head[u]

    C = tuple(range(r, r + m))
    X = tuple(range(r + m, r + 2 * m))

    # C is a clique.
    for u, v in itertools.combinations(C, 2):
        add_edge(adj, u, v)

    # R is complete to C and anticomplete to X.
    for u in range(r):
        for c in C:
            add_edge(adj, u, c)

    # C-X adjacency.
    if kind == "thin":
        for i in range(m):
            add_edge(adj, C[i], X[i])
    elif kind == "thick":
        for i in range(m):
            for j in range(m):
                if i != j:
                    add_edge(adj, C[i], X[j])
    else:
        raise ValueError(f"unknown spider kind: {kind}")

    return tuple(adj), C, X


def duplicate_vertex(graph: Adj, v: int, twin_type: str) -> Adj:
    """Replace/duplicate v by adding a true or false twin of v."""
    n = len(graph)
    adj = list(graph) + [0]
    old_neighbors = graph[v]

    # New vertex has every old neighbor of v.
    bitset = old_neighbors
    while bitset:
        lsb = bitset & -bitset
        u = lsb.bit_length() - 1
        add_edge(adj, n, u)
        bitset ^= lsb

    if twin_type == "true":
        add_edge(adj, n, v)
    elif twin_type != "false":
        raise ValueError(f"unknown twin type: {twin_type}")

    return tuple(adj)


def quasi_spider(
    head: Adj,
    m: int,
    kind: str,
    duplicated_side: str,
    twin_type: str,
    s: int = 0,
) -> Adj:
    graph, C, X = spider(head, m, kind)
    v = C[s] if duplicated_side == "C" else X[s]
    return duplicate_vertex(graph, v, twin_type)


def alpha_vector(graph: Adj) -> List[int]:
    """Compute alpha_k(G) exactly for 1 <= k <= |V(G)|.

    A set S is k-independent when Delta(G[S]) <= k-1.
    Exhaustive enumeration is feasible here because the default tests
    only reach order 13.
    """
    n = len(graph)
    if n == 0:
        return [0]

    # best_exact[d] = largest |S| encountered with Delta(G[S]) == d.
    # Use index 0 for d=-1 (the empty set), d+1 otherwise.
    best_exact = [0] * (n + 1)

    for mask in range(1 << n):
        size = mask.bit_count()
        if mask == 0:
            dmax = -1
        else:
            dmax = 0
            bits = mask
            while bits:
                lsb = bits & -bits
                v = lsb.bit_length() - 1
                deg = (graph[v] & mask).bit_count()
                if deg > dmax:
                    dmax = deg
                bits ^= lsb
        idx = dmax + 1
        if size > best_exact[idx]:
            best_exact[idx] = size

    # alpha_k is the maximum size among subsets with dmax <= k-1.
    alpha = [0] * (n + 1)
    running = 0
    for k in range(1, n + 1):
        # admissible dmax values are -1,0,...,k-1 -> indices 0,...,k
        running = max(running, best_exact[k])
        # best_exact[0] is 0, so it never changes running.
        # We also need to accumulate all previous exact values.
        if k == 1:
            running = max(best_exact[0], best_exact[1])
        else:
            running = max(alpha[k - 1], best_exact[k])
        alpha[k] = running
    return alpha


def alpha_head(head: Adj, max_k: int) -> List[int]:
    """alpha_k(G[R]) for all k needed by the formulas."""
    r = len(head)
    values = [0] * (max_k + 1)
    if r == 0:
        return values

    exact = alpha_vector(head)
    for k in range(1, max_k + 1):
        values[k] = exact[k] if k <= r else r
    return values


def formula_value(
    structure: Structure,
    k: int,
    m: int,
    r: int,
    a: Sequence[int],
) -> int:
    """Closed formulas from Table 1 of the manuscript."""
    ak = a[k]

    if structure.family == "spider":
        if structure.kind == "thin":
            return m + max(ak, min(k - 1, m + r))
        if structure.kind == "thick":
            if m <= k and r < k - m + 1:
                return min(k + 1, 2 * m + r)
            return m + ak
        raise ValueError("unknown spider kind")

    kind = structure.kind
    side = structure.duplicated_side
    twin = structure.twin_type

    # Thin (R,C,X^f)
    if kind == "thin" and side == "X" and twin == "false":
        if k == 1:
            return m + 1 + a[1]
        return m + 1 + max(
            ak,
            min(k - 1, r + m - 1),
            min(k - 2, r + m),
        )

    # Thin (R,C,X^t)
    if kind == "thin" and side == "X" and twin == "true":
        if k == 1:
            return m + a[1]
        return m + 1 + max(
            ak,
            min(k - 1, r + m - 1),
            min(k - 2, r + m),
        )

    # Thin (R,C^t,X)
    if kind == "thin" and side == "C" and twin == "true":
        return m + max(ak, min(k - 1, r + m + 1))

    # Thin (R,C^f,X)
    if kind == "thin" and side == "C" and twin == "false":
        if k <= 2:
            return m + max(ak, 1)
        return m + max(
            ak,
            min(k - 1, r + m + 1),
            min(k, r + 2),
        )

    # Thick (R,C,X^f)
    if kind == "thick" and side == "X" and twin == "false":
        if k < m:
            return m + 1 + ak
        return m + 1 + max(
            ak,
            min(k - m, r + m),
            min(k - m + 1, r + 1),
        )

    # Thick (R,C,X^t)
    if kind == "thick" and side == "X" and twin == "true":
        if k == 1:
            return m + a[1]
        if 2 <= k < m:
            return m + 1 + ak
        return m + 1 + max(
            ak,
            min(k - m, r + m),
            min(k - m + 1, r + 1),
        )

    # Thick (R,C^t,X)
    if kind == "thick" and side == "C" and twin == "true":
        if k < m:
            return m + ak
        return m + max(ak, min(k - m + 1, r + m + 1))

    # Thick (R,C^f,X)
    if kind == "thick" and side == "C" and twin == "false":
        if k <= 2:
            return max(m + ak, 3)
        akm2 = a[k - 2]
        if 3 <= k < m:
            return max(
                m + ak,
                3 + min(k - 1, m - 1 + akm2),
            )
        return max(
            m + ak,
            m + min(k - m + 1, r + m + 1),
            m + 2 + min(r, k - m),
            3 + min(k - 1, m - 1 + akm2),
        )

    raise ValueError(f"unrecognized structure: {structure}")


def structure_types() -> Tuple[Structure, ...]:
    structures: List[Structure] = [
        Structure("spider", "thin"),
        Structure("spider", "thick"),
    ]
    for kind in ("thin", "thick"):
        for side in ("X", "C"):
            for twin in ("false", "true"):
                structures.append(Structure("quasi", kind, side, twin))
    return tuple(structures)


def run_tests(max_m: int, max_r: int, verbose: bool = False) -> Tuple[int, int]:
    structures = structure_types()
    instances = 0
    comparisons = 0

    for m in range(2, max_m + 1):
        for r in range(0, max_r + 1):
            for head_index, head in enumerate(all_graphs(r)):
                max_needed = r + 2 * m + 1
                a = alpha_head(head, max_needed)

                for structure in structures:
                    if structure.family == "spider":
                        graph, _, _ = spider(head, m, structure.kind)
                    else:
                        graph = quasi_spider(
                            head,
                            m,
                            structure.kind,
                            structure.duplicated_side or "",
                            structure.twin_type or "",
                        )

                    instances += 1
                    exact = alpha_vector(graph)
                    n = len(graph)

                    for k in range(1, n + 1):
                        predicted = formula_value(structure, k, m, r, a)
                        comparisons += 1
                        if predicted != exact[k]:
                            print("MISMATCH", file=sys.stderr)
                            print(f"structure: {structure.label}", file=sys.stderr)
                            print(f"m={m}, r={r}, head_index={head_index}, k={k}", file=sys.stderr)
                            print(f"formula={predicted}, exact={exact[k]}", file=sys.stderr)
                            print(f"head adjacency masks={head}", file=sys.stderr)
                            raise SystemExit(1)

                if verbose:
                    print(f"checked m={m}, r={r}, head #{head_index}")

    return instances, comparisons


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exhaustively verify the spider/quasi-spider formulas in Table 1."
    )
    parser.add_argument(
        "--max-m",
        type=int,
        default=4,
        help="largest spider weight m to test (default: 4)",
    )
    parser.add_argument(
        "--max-r",
        type=int,
        default=4,
        help="largest head order |R| to test (default: 4)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print progress for each generated head",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_m < 2:
        raise SystemExit("--max-m must be at least 2")
    if args.max_r < 0:
        raise SystemExit("--max-r must be nonnegative")

    instances, comparisons = run_tests(args.max_m, args.max_r, args.verbose)
    print("All tests passed.")
    print(f"Structural instances: {instances}")
    print(f"Comparisons: {comparisons}")


if __name__ == "__main__":
    main()
