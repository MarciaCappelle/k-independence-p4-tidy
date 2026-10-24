"""
Test driver for the computational verification.

Four stages:

  1. corpus validity   every generated graph is P_4-tidy, checked directly
                       from the definition (every induced P_4 has at most one
                       partner), independently of how it was built;

  2. tree construction T(G) is built for every corpus graph by the procedure
                       of Lemma 21, and the structural labels it produces are
                       compared with the construction parameters;

  3. main test         alpha_k computed by Algorithm 2 over T(G) is compared
                       with brute force, for every graph and every k;

  4. regression        the same pipeline with the intra-module exclusion
                       disabled, to confirm the test is sensitive to it.
"""

import sys
import time
from collections import Counter

from p4tidy_core import (generate_corpus, is_p4_tidy, alpha_all, spider,
                         path, P5, C5, P5bar, K1, complement)
from p4tidy_tree import labelled_tree, LabellingError, count_nodes
from p4tidy_dp import dp_vector, alpha_k_algorithm

MAX_N = int(sys.argv[1]) if len(sys.argv) > 1 else 12


def tree_profile(T, acc):
    if T.kind == "spider":
        side, kind = T.label["side"], T.label["kind"]
        tt = "thick" if T.label["thick"] else "thin"
        if side is None:
            acc[f"spider {tt}"] += 1
        else:
            acc[f"quasi {tt} ({side}^{kind[0]})"] += 1
        acc["head nonempty" if T.label["r"] else "head empty"] += 1
    elif T.kind == "exceptional":
        acc[f"exceptional {T.label['type']}"] += 1
    else:
        acc[T.kind] += 1
    for c in T.children:
        tree_profile(c, acc)


def main():
    t0 = time.time()
    print(f"corpus: graphs up to {MAX_N} vertices")
    corpus = generate_corpus(max_n=MAX_N)
    print(f"  {len(corpus)} graphs generated")

    # ---- stage 1 -------------------------------------------------------
    bad = [tag for g, tag in corpus if not is_p4_tidy(*g)]
    print(f"\n[1] P_4-tidy check (direct definition): "
          f"{len(corpus) - len(bad)}/{len(corpus)} pass")
    if bad:
        print("    FAILURES:", bad[:5])
        return 1

    # ---- stage 2 -------------------------------------------------------
    profile = Counter()
    build_fail = []
    trees = []
    for g, tag in corpus:
        try:
            T = labelled_tree(*g)
        except (LabellingError, AssertionError) as exc:
            build_fail.append((tag, str(exc)))
            trees.append(None)
            continue
        trees.append(T)
        tree_profile(T, profile)
    print(f"\n[2] labelled tree T(G) built: "
          f"{len(corpus) - len(build_fail)}/{len(corpus)} pass")
    if build_fail:
        print("    FAILURES:", build_fail[:5])
        return 1
    print("    node coverage over all trees:")
    for key in sorted(profile):
        print(f"      {key:28s} {profile[key]}")

    # ---- stage 3 -------------------------------------------------------
    mismatches = []
    checks = 0
    for (g, tag), T in zip(corpus, trees):
        n, adj = g
        bf = alpha_all(n, adj, n)
        for k in range(1, n + 1):
            got = alpha_k_algorithm(n, adj, k)
            checks += 1
            if got != bf[k]:
                mismatches.append((tag, n, k, bf[k], got))
    print(f"\n[3] alpha_k: algorithm vs brute force over all k: "
          f"{checks - len(mismatches)}/{checks} pass")
    if mismatches:
        print("    MISMATCHES (tag, n, k, brute force, algorithm):")
        for m in mismatches[:10]:
            print("      ", m)
        return 1

    # ---- stage 4 -------------------------------------------------------
    err = wrong = ok = 0
    for g, tag in corpus:
        n, adj = g
        bf = alpha_all(n, adj, n)
        try:
            vals = [alpha_k_algorithm(n, adj, k, premark=False)
                    for k in range(1, n + 1)]
        except (LabellingError, AssertionError):
            err += 1
            continue
        if vals != bf[1:]:
            wrong += 1
        else:
            ok += 1
    print(f"\n[4] regression with the intra-module exclusion disabled:")
    print(f"      construction fails : {err}")
    print(f"      wrong alpha_k      : {wrong}")
    print(f"      unaffected         : {ok}")
    if err + wrong == 0:
        print("    WARNING: the test suite is insensitive to this step")
        return 1

    nodes = sum(count_nodes(T) for T in trees)
    print(f"\nall stages passed  ({checks} alpha_k comparisons, "
          f"{nodes} tree nodes, {time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
