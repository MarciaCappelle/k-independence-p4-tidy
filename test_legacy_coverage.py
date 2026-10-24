"""Exhaustive labelled-head regression for Table 1 (original verification scope).

Run separately: python3 test_legacy_coverage.py
This deliberately tests formulas even for heads that are not P4-tidy.
"""
from itertools import combinations
from p4tidy_core import spider, alpha_all, add_edge
from p4tidy_dp import spider_formula


def labelled_heads(r):
    edges = list(combinations(range(r), 2))
    for mask in range(1 << len(edges)):
        adj = [0] * r
        for i, (u, v) in enumerate(edges):
            if (mask >> i) & 1:
                add_edge(adj, u, v)
        yield r, adj


def main():
    instances = comparisons = 0
    for m in (2, 3, 4):
        for r in range(5):
            for head_n, head_adj in labelled_heads(r):
                head = (head_n, head_adj) if r else None
                for thick in (False, True):
                    for side, kind in [(None, None), ("X", "false"), ("X", "true"),
                                       ("C", "false"), ("C", "true")]:
                        n, adj = spider(m, thick, head, side, kind)
                        exact = alpha_all(n, adj, n)
                        a = ([0] + [0] * n) if r == 0 else alpha_all(r, head_adj, n)
                        label = {"m": m, "r": r, "thick": thick,
                                 "side": side, "kind": kind}
                        instances += 1
                        for k in range(1, n + 1):
                            got = spider_formula(label, k, a)
                            comparisons += 1
                            if got != exact[k]:
                                print(f"FAIL: m={m} r={r} thick={thick} side={side} "
                                      f"kind={kind} head={head_adj} k={k} "
                                      f"formula={got} exact={exact[k]}")
                                return 1
    print("All legacy-scope formula tests passed.")
    print(f"Structural instances: {instances}")
    print(f"Comparisons: {comparisons}")
    if (instances, comparisons) != (2280, 24054):
        print("WARNING: unexpected historical coverage counts")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
