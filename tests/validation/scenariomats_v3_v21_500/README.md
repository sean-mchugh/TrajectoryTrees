# Preserved 500-tree similarity equality audit

This directory preserves the fixed 500-tree differential audit originally used
to compare ScenarioMats V3 against the plain-R ScenarioMats V2.1 numerical
authority. Within TrajectoryTrees, the former V3 candidate is exposed as the
version-neutral public function `TrajectoryTrees::trajectory_similarity()`.

## Corpus

`selected_inventory.csv` is the frozen 500-row selection manifest. The original
source paths remain in the manifest as provenance. Self-contained copies are
stored under `trees/` with the same RDS basenames.

The fixed corpus contains:

- 420 adversarial and 80 generated stochastic character maps;
- 322 ultrametric and 178 non-ultrametric trees;
- trees with 3, 4, 5, 6, 8, 9, or 10 tips;
- between 1 and 7 observed states.

## What the audit does

For every selected simmap, the runner:

1. derives the construction tolerance from the smallest positive mapped
   duration, capped at `1e-10`;
2. constructs one `trajectory_obj` with
   `TrajectoryTrees::make_trajectory_objects()` and retains scenario matrices
   and path maps;
3. passes that exact object independently to the test-only V2.1 oracle and
   `TrajectoryTrees::trajectory_similarity()`;
4. compares every candidate matrix and reported summary against V2.1 for depths
   0-3, both `bhattacharyya` and `minimum`, and both `complete` and
   `matrices_and_summaries` return modes at tolerance `1e-10`;
5. checks matrix dimensions, labels, finiteness, non-negativity, symmetry,
   zero diagonals, and pairwise closure;
6. independently reconstructs each summary's totals, tree-wide proportions,
   means, total across matrices, and—when non-ultrametric weighting is
   active—available-similarity total from the corresponding returned matrices;
7. checks tree-wide closure, correct ultrametric/non-ultrametric weighting, and
   complete-versus-compact consistency.

The package-only `summaries$similarity_vectors` extension is tested elsewhere
and is removed only from the V2.1 equality surface because V2.1 predates that
extension.

## Run it

Start in the TrajectoryTrees package root:

```text
/Users/seanmchugh/Projects/TrajectoryTrees
```

Ensure TrajectoryTrees 0.2.0 or newer is installed, then run:

```bash
Rscript tests/validation/scenariomats_v3_v21_500/run_audit.R
```

The script takes no trailing arguments. A successful full run ends with:

```text
TrajectoryTrees similarity equality passed for 500 histories.
```

## Results

Every invocation creates a uniquely named directory under `results/`. Each
result directory contains:

- `audit_metadata.csv`: status, corpus size, package version, calculation
  surfaces, numerical tolerance, and working directory;
- `discrepancies.csv`: every V2.1/package-native structural or numerical
  discrepancy; it is header-only on success;
- `internal_checks.csv`: the invariant result for every history, metric, return
  mode, and matrix family.

The runner stops with a nonzero status when any discrepancy or internal check
fails and prints the result directory containing the evidence.
