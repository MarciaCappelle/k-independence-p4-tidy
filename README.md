# Verification of the k-independence formulas for P4-tidy graphs

This repository contains the exhaustive verification script associated with the manuscript **"k-independence in P4-tidy graphs"**.

The script independently checks the closed formulas for the k-independence number of the spider and quasi-spider structures summarized in Table 1 of the manuscript.

## What is verified

For every generated instance, the script:

1. constructs the relevant spider or quasi-spider;
2. computes the exact value of `alpha_k(G)` by exhaustive enumeration of all vertex subsets;
3. evaluates the corresponding closed formula from Table 1; and
4. compares the two values for every `1 <= k <= |V(G)|`.

The default run considers:

- spider weight `m` in `{2, 3, 4}`;
- head order `0 <= |R| <= 4`;
- **every labelled simple graph** on the head `R`;
- the two spider types: thin and thick;
- the eight quasi-spider types obtained by duplicating a vertex of `C` or `X` as true or false twins, in both the thin and thick cases;
- every value `1 <= k <= |V(G)|`.

The largest generated graphs have order 13.

With the default parameters, the script checks **2,280 structural instances** and performs **24,054 formula-versus-exact comparisons**.

## Requirements

- Python 3.8 or later
- No third-party packages

## Running the verification

```bash
python verify_formulas.py
```

Expected output:

```text
All tests passed.
Structural instances: 2280
Comparisons: 24054
```

A nonzero exit status is returned if a discrepancy is found. In that case, the script prints the structural type, `m`, `|R|`, the index of the labelled head graph, `k`, the value predicted by the formula, the exact value, and the adjacency masks of the head.

## Changing the verification range

The default limits reproduce the verification reported in the revised manuscript. They can be changed from the command line. For example:

```bash
python verify_formulas.py --max-m 4 --max-r 5
```

The number of generated heads grows rapidly: there are `2^(r choose 2)` labelled simple graphs on `r` vertices. Consequently, increasing `--max-r` is considerably more expensive than increasing `--max-m`.

Use `--verbose` to print progress information:

```bash
python verify_formulas.py --verbose
```

## Scope of the computation

This computation is an **independent consistency check** of the formulas and their boundary cases. It is not used as a substitute for the mathematical proofs in the manuscript.

The head `R` is not restricted to be P4-tidy in this verification. Every simple graph on at most four vertices is used. This is stronger for the local formulas, which depend on the head through its order and relevant k-independence numbers.

## Reproducibility

The script generates all tested instances internally; no external graph-instance files or datasets are required.

For archival citation, the intended workflow is to create a tagged GitHub release and archive that release with Zenodo. The Zenodo DOI can then be cited in the revised manuscript.
