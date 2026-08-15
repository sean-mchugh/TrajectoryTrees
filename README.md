# TrajectoryTrees

TrajectoryTrees is an R package for turning stochastic character maps into
three complementary evolutionary representations:

- the mapped phylogeny;
- a scenario tree that separates distinct lineage histories;
- a transition tree that summarizes realized character-state sequences.

It also constructs summary-statistic and through-time surfaces, aggregates
ensembles into posterior distributions, writes large jobs in resumable batches,
and provides focused plotting functions.

## Installation

From the parent directory:

```r
install.packages("TrajectoryTrees", repos = NULL, type = "source")
```

Or during development:

```r
pkgload::load_all("TrajectoryTrees")
```

Package installation compiles the native implementation once using C++17.
Ordinary use does not invoke a compiler or call `sourceCpp()`.

## Constructing trajectory objects

```r
library(TrajectoryTrees)

simmap <- readRDS(system.file(
  "extdata", "generated_000001.rds",
  package = "TrajectoryTrees"
))

trajectory <- make_trajectory_objects(simmap)
trajectory
```

The same function accepts a named list or a named list of lists and preserves
that shape:

```r
posterior <- make_trajectory_objects(list(map_1 = simmap, map_2 = simmap))

groups <- make_trajectory_objects(list(
  empirical = list(map_1 = simmap),
  simulated = list(map_1 = simmap)
))
```

Keep a compact result when only selected surfaces are needed:

```r
compact <- select_trajectory_fields(
  trajectory,
  include = c("phylo", "scenario", "trans", "ss")
)
```

## Summaries and distributions

```r
state_statistics <- trajectory_summary(trajectory, "by_state")
lineages <- trajectory_through_time(trajectory, "ltt", "tot")

distribution <- make_trajectory_distributions(
  posterior,
  include = c("P", "T", "SS", "TT")
)
```

`compare_trajectory_distributions()` recursively compares matching numeric
surfaces and reports the maximum absolute difference for each path.

## Plotting

```r
plot_trajectory_tree(trajectory, tree = "scenario")

plot_trajectory_through_time(
  trajectory,
  metric = "ltt",
  component = "state"
)

plot_trajectory_distribution(distribution, surface = "SS")
```

`plot_trajectory_tanglegram()` compares two tree representations. It uses
`phytools` when available and accepts an explicit two-column tip association
matrix for specialized labels.

## Resumable batches

```r
manifest <- make_trajectory_batches(
  posterior_maps,
  output_dir = "trajectory_batches",
  batch_size = "auto"
)

posterior <- load_trajectory_batches(manifest)
```

Completed batch files are never silently overwritten. An interrupted compatible
run can continue with `resume = TRUE`.

## Bundled validation and Anolis trees

The source package includes the exact release datasets:

- 626 validation trees: 471 ultrametric and 155 non-ultrametric;
- 500 empirical Anolis stochastic maps;
- 1,000 simulated Anolis stochastic maps.

Load or inspect a filtered validation subset before reading any trees:

```r
inventory <- load_validation_trees(
  min_tips = 10,
  max_tips = 100,
  exclude_label = "^generated_",
  load = FALSE
)

trees <- load_validation_trees(
  population = "nonultrametric",
  exclude_tips = c(2, 3),
  limit = 25
)

empirical <- load_anolis_trees("empirical")
simulated <- load_anolis_trees("simulated", limit = 50)
```

## Known-true references

Known-true cases receive a separate validation gate. The repository contains
129 registered saved references:

| Scope | Cases |
|---|---:|
| PST | 20 |
| PST + SS | 39 |
| PST + TT + SS | 21 |
| Original metrics | 49 |

The complete inventory is
`tools/validation/reference/true/manifest_complete.csv`. The former 91-row
manifest is retained for provenance; none of its saved values were rewritten.

Run the known-true gate from the package directory:

```sh
Rscript tools/validation/run_validation.R \
  --scopes=known-true \
  --output=tools/validation/output/known_true
```

## Validation and the historical reference implementation

The production implementation is private and all exported functions use
version-neutral names. V34 is retained only under
`tools/validation/reference/v34/` as a historical comparison oracle. It is not
compiled into the package, loaded by normal use, exported, or called by any
public function.

Optional hidden parity:

```sh
Rscript tools/validation/run_validation.R \
  --scopes=corpus,anolis-empirical,anolis-simulated \
  --compare-v34=true \
  --output=tools/validation/output/full_parity
```

See `tools/validation/README.md` for every filter. Compact historical release
summaries are retained under `tools/validation/historical/`.

## Repository layout

```text
R/                    public API and private R implementation
src/                  package-native C++ implementation
man/                  generated function documentation
tests/testthat/       fast package tests
vignettes/            core and Anolis workflows
inst/extdata/         examples and release tree datasets
tools/validation/     full gates, known truths, and hidden V34
```

## License

Copyright 2026 Sean McHugh. All rights reserved. See `LICENSE`.
