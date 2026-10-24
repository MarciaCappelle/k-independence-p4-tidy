# Computational verification for k-independence in P4-tidy graphs

This repository contains the computational verification code associated with the manuscript **"k-independence in P4-tidy graphs"**.

The code verifies the construction of the labelled decomposition tree `T(G)` (Lemma 21), the closed formulas for spiders and quasi-spiders (Table 1), and the dynamic programming algorithm (Algorithms 1 and 2).

The results produced by the algorithm are compared with exact values obtained by exhaustive enumeration of vertex subsets.

## Archived version

The verification code is archived on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.23619639

The repository is being updated to include the complete verification suite. A new Zenodo version will be published after the updated GitHub release.

## Files

- `p4tidy_core.py`: graph representation, graph constructors, corpus generation, direct P4-tidy recognition, and brute-force computation of `alpha_k`.
- `p4tidy_tree.py`: modular decomposition, construction of the labelled tree `T(G)`, and computation of representative-graph degrees.
- `p4tidy_dp.py`: closed formulas from Table 1 and the dynamic programming algorithms.
- `test_all.py`: main verification program.
- `test_log.txt`: output of the complete verification run.
- `test_legacy_coverage.py`: exhaustive verification of the Table 1 formulas for all labelled head graphs of order at most four.

## Requirements

- Python 3.8 or later
- No third-party packages

## Running the verification

Run the following command:

```bash
python3 test_all.py 13
```

The argument `13` specifies the maximum order of the generated graphs.

For a smaller verification run:

```bash
python3 test_all.py 11
```

The program reports the number of successful tests, the structural types encountered, and any discrepancies.

## What is verified

The verification consists of four stages.

### 1. Corpus validity

Graphs are generated using the structural characterization of P4-tidy graphs, including disjoint unions, joins, spiders, quasi-spiders, and exceptional graphs.

Every generated graph is independently checked against the definition of P4-tidy graphs.

### 2. Construction of the labelled decomposition tree

For every generated graph, the program constructs its modular decomposition and the labelled decomposition tree `T(G)`.

The verification covers:

- thin and thick spiders;
- all eight types of quasi-spiders;
- empty and nonempty heads;
- nested decomposition structures;
- the exceptional graphs `P5`, `P5bar`, and `C5`;
- modules consisting of true or false twins.

The representative-graph degrees are computed by scanning adjacency lists and counting distinct adjacent children.

### 3. Verification of the k-independence algorithm

For every generated graph `G` and every integer `1 <= k <= |V(G)|`, the program:

1. computes the exact value of `alpha_k(G)` by exhaustive enumeration of vertex subsets;
2. constructs the labelled decomposition tree `T(G)`;
3. evaluates the dynamic programming algorithm; and
4. compares the algorithmic result with the exact value.

### 4. Regression test for Lemma 21

The verification also tests the exclusion of edges internal to a module when computing degrees in the representative graph.

If a representative belongs to a two-vertex series module, its adjacent twin must not be counted as a neighboring child.

The program repeats the construction with this exclusion disabled to verify that the test suite detects the resulting errors.

## Computational results

For the corpus containing graphs of order at most 13, the recorded verification produced:

- **7,421 generated graphs**;
- **83,904 comparisons** between the algorithm and brute force;
- **42,883 decomposition-tree nodes**;
- **no discrepancies** in the computed k-independence numbers.

In the regression test, disabling the exclusion of intra-module edges caused construction failures in **4,280 graphs**.

The complete recorded output is available in `test_log.txt`.
## Additional exhaustive verification of Table 1

The script `test_legacy_coverage.py` independently verifies the closed formulas for thin and thick spiders and all eight quasi-spider types.

It considers spider weights m = 2, 3, 4 and every labelled simple graph on a head R of order 0 to 4.

The verification covers 2,280 structural instances and 24,054 comparisons with brute-force values, with no discrepancies.

Run it using:

```bash
python3 test_legacy_coverage.py
```

Unlike the main corpus, this verification allows arbitrary head graphs, including graphs that are not P4-tidy.

## Scope and limitations

The verification provides an independent computational consistency check of the structural formulas, decomposition-tree construction, and dynamic programming algorithm.

The generated corpus includes a broad range of structural configurations but is not an exhaustive enumeration of all P4-tidy graphs up to 13 vertices.

The verification implementation uses exhaustive module enumeration to construct the modular decomposition. It is intended for testing correctness, not for achieving the linear-time modular decomposition bound assumed in the theoretical complexity analysis.

The computational results do not replace the mathematical proofs presented in the manuscript.

## Reproducibility

All graph instances are generated internally. No external datasets or third-party packages are required.

The test suite can be reproduced by running `test_all.py` with the desired maximum graph order.

The corresponding GitHub release will be archived on Zenodo to preserve the version associated with the manuscript.
