# TrajectoryTrees test and validation map

Unless a section explicitly names the `Scenario_tree_fast_codex` root, run its
commands from the TrajectoryTrees package root:

```text
/Users/seanmchugh/Projects/TrajectoryTrees
```

## Test layers are not interchangeable

There are three distinct validation layers:

1. the authoritative trajectory-calculation contracts and large reference
   campaigns;
2. lightweight package-interface tests under `tests/testthat/`; and
3. ScenarioMats matrix and summary validation performed only after a
   `trajectory_obj` has been constructed.

The ordinary `testthat` suite is **not** the numerical correctness authority
for `trajectory_obj`. It checks the installed public interface and object
handling. Trajectory-calculation correctness is established in three steps:

1. independently recover the scenario matrices from the original simmap;
2. independently recover the scenario-tree and transition-tree quantities from
   those scenario matrices; and
3. compare complete PST, SS, and TT results with stored, by-hand-confirmed
   reference objects.

## 1. Authoritative `trajectory_obj` correctness validation

The full pre-package V35 correctness framework remains at:

```text
/Users/seanmchugh/Projects/Scenario_tree_fast_codex/testing/correctness/final_tests/
```

These tests construct real trajectory objects and verify their calculated
scenario, transition, SS, and TT contents. They do not merely check class names,
field names, or dimensions.

### A small example: movie, spreadsheet, and summaries

Imagine that a simmap is a movie of colored branches growing through time.
State `A` is one color and state `B` is another. The movie contains three tips:

- `tip_1` stays in state `A` for its entire history;
- `tip_2` changes from `A` to `B` at time 3; and
- `tip_3` changes from `A` to `B` at time 4.

All three tips finish at time 5. A simplified state-only view of the movie is:

| Tip | Time 0 | Time 2 | Time 3 | Time 4 | Time 5 |
|---|---|---|---|---|---|
| `tip_1` | A | A | A | A | A |
| `tip_2` | A | A | B | B | B |
| `tip_3` | A | A | A | B | B |

The real `scenario_mats` object contains this state grid plus several parallel
grids. At every populated cell it records:

- which original phylogeny edge the tip occupies;
- which mapped step on that phylogeny edge it occupies;
- which scenario-tree edge and mapped step it occupies;
- its complete phylogeny path; and
- its complete scenario path.

It is useful to think of the simmap as the original movie and `scenario_mats`
as a frame-by-frame spreadsheet made from that movie. The scenario tree and
transition tree are then two different summaries made from the spreadsheet:

- the **scenario tree** keeps distinct lineage histories separate; and
- the **transition tree** groups histories by their realized character-state
  sequence.

In the example, `tip_2` and `tip_3` both realize the sequence `A -> B`, so that
sequence is represented together in the transition summary. Their changes
happen on different lineages and at different times, however, so their detailed
scenario histories remain distinguishable.

The validation deliberately checks this chain in two directions:

```text
original simmap  --->  scenario_mats  --->  scenario and transition trees
       step 1                 step 2
```

Step 1 asks whether the spreadsheet faithfully copied the movie. Step 2 treats
the spreadsheet as evidence and asks whether the two summaries were calculated
correctly from it. A mistake cannot be hidden merely because two derived
objects happen to agree with each other.

### A. Does `scenario_mats` reproduce the original simmap?

For this step, the test does not trust the coordinates already stored in
`scenario_mats`. It independently walks from the root through every edge in the
simmap. On each edge it adds the mapped durations to recover the absolute start
and end time of every state step.

For example, suppose the original map for `tip_2` says:

```text
A from time 2 through time 3
B from time 3 through time 5
```

The first matrix cell that says `B` for `tip_2` must therefore occur at time 3,
point to the correct phylogeny edge, and point to the `B` step on that edge. The
test follows those stored IDs back into the simmap and asks, in effect:

> If I use this cell's edge number and step number to look at the original
> movie, do I recover exactly the state, path, and time written in the cell?

It repeats that question for every populated cell, for both phylogeny and
scenario coordinates. The test fails if, for example:

- `B` appears one time column too early or too late;
- an edge ID names the wrong branch;
- a step ID names the wrong piece of a branch;
- a stored state differs from the simmap state; or
- a stored path differs from the path obtained from the simmap.

| Test family | Independent check |
|---|---|
| `scenario_mats_edge_step_times` | Reconstructs the absolute start and end time of every simmap step from the tree and its mapped durations. Every scenario-matrix edge/step coordinate must occur inside the corresponding simmap interval and begin at the correct time. |
| `scenario_mats_edge_steps` | Uses each matrix edge ID and step ID to look the value up directly in the original scenario-tree and phylogeny simmaps. The recovered state and path must exactly equal the state and path stored in that matrix cell. |

These checks answer whether the rectangular matrices are an accurate rendering
of the supplied mapped history.

### B. Do the scenario and transition trees agree with `scenario_mats`?

For this step, the test temporarily ignores the sizes and classifications
already stored on the scenario and transition trees. It recalculates them from
the matrix rows.

#### Rebuilding the scenario-tree information

The test reads each tip's matrix row from left to right. Tips belong to the same
scenario history only while their complete histories agree. In the example:

- all three tips initially share history `A`;
- `tip_2` leaves that history when it changes to `B` at time 3; and
- `tip_3` leaves it when it changes to `B` at time 4.

For every scenario edge and path segment, the test finds all matrix cells that
claim membership in it. It then independently checks:

1. all included rows really have the same state history through that point;
2. no excluded row has that same history and was accidentally left out;
3. the stored path is the path spelled out by those matrix states; and
4. the stored lineage size equals the number of distinct phylogeny lineages
   represented by those rows.

For a non-ultrametric tree, tips do not all finish at the same time. The test
therefore repeats the lineage count through time, checking the complete stored
size trajectory instead of checking only one final count.

#### Rebuilding the transition-tree information

The transition tree answers a different question: how many lineages and
scenario histories realize each state-transition path? In the example,
`tip_2` and `tip_3` both contribute to the `A -> B` path even though their
transitions occur at different times.

For each transition-tree edge and path step, the test:

1. finds the phylogeny tips represented below that transition edge;
2. uses their matrix rows to locate the relevant path;
3. finds the terminal contributing phylogeny edges and scenario edges;
4. counts those contributors independently; and
5. compares the counts with the lineage-size and scenario-size maps stored on
   the transition tree.

Finally, after the scenario paths themselves have been checked against the
matrices, a separate topology check ignores the stored `edge_type` labels. It
looks at each scenario edge's parent:

- the edge with no parent is `Root`;
- a child whose first path equals its parent's final path is `Stay`; and
- a child whose first path differs from its parent's final path is `Leave`.

Those independently derived labels must equal the stored labels. The test then
counts the derived `Leave` edges entering each destination path and compares
those counts with the `lindiff` values in the SS by-path table.

The key point is that the test does not say, "the scenario tree and transition
tree look reasonable." It derives the expected answers again from the matrix
evidence and compares the two calculations.

| Test family | Independent check |
|---|---|
| `scenario_edge_state_paths` | Reconstructs the state histories belonging to each scenario edge directly from the matrix rows and checks the stored scenario-edge state paths. |
| `scenario_edge_paths` | Reconstructs scenario membership, path runs, and final lineage ownership from the matrices and checks the scenario tree's stored paths and sizes. |
| `scenario_edge_lineage_size_trajectories` | For non-ultrametric trees, reconstructs the lineage count through time for each scenario edge from matrix ownership and checks every stored size trajectory. |
| `transition_edge_path_sizes` | Finds the matrix rows represented below each transition-tree edge and independently counts their contributing phylogeny and scenario edges. Those counts must equal the transition tree's stored lineage and scenario sizes. |
| `scenario_edge_stay_leave_types` | Reconstructs whether lineages stay on or leave each scenario path and checks the stored event classifications. |
| `scenario_leave_counts_by_path` | Counts departures from each path directly from the matrix histories and checks the stored leave counts. |

The ultrametric suite contains the two matrix-versus-simmap checks and five
scenario/transition checks. The non-ultrametric suite adds the explicit
lineage-size-through-time check, giving eight test families.

The implementations are:

| File | Contents |
|---|---|
| `src/ultrametric/pst_reference_tests_ultrametric.R` | Registers and runs the seven ultrametric test families. |
| `src/ultrametric/ultra_matrix_coordinates.R` | Matrix-versus-simmap time, edge, step, state, and path checks. |
| `src/ultrametric/ultra_scenario_edges.R` and `ultra_transition_sizes.R` | Matrix-derived scenario-tree and transition-tree checks. |
| `src/nonultrametric/pst_reference_tests_nonultrametric.R` | Registers and runs the eight non-ultrametric test families. |
| `src/nonultrametric/nonultra_matrix_coordinates.R`, `nonultra_scenario_edges.R`, and `nonultra_transition_sizes.R` | Endpoint-aware versions of the matrix, scenario-tree, and transition-tree checks. |
| `release/run_v35_release_barrier.R` | Applies the appropriate seven- or eight-family suite to all 626 adversarial/generated trees: 471 ultrametric and 155 non-ultrametric. It also requires successful construction, unchanged input simmaps, and correct ultrametric/non-ultrametric identification. |

The 626-tree collection is a broad input corpus for these independent logic
tests. It is not a collection of 626 stored expected trajectory objects.

Run the direct 626-tree V35 logic-test barrier from the
`Scenario_tree_fast_codex` root:

```bash
Rscript testing/correctness/final_tests/release/run_v35_release_barrier.R \
  --output=testing/correctness/final_tests/evidence/v35_main_release_candidate
```

### C. Does the result match by-hand-confirmed PST/SS/TT references?

The TrajectoryTrees source checkout contains 129 stored known-true reference
objects under `tools/validation/reference/true/`. The manifest specifies which
parts of each object were confirmed and must be compared:

| Confirmed reference scope | Number of objects |
|---|---:|
| PST structure | 20 |
| PST and SS | 39 |
| PST, TT, and SS | 21 |
| Original SS metrics | 49 |
| **Total** | **129** |

For a PST reference, the runner compares `phylo`, `scenario`, `trans`, and
`root_policy`. It additionally compares `ss` and `tt` when those fields are
marked as confirmed in `manifest_complete.csv`.

The TrajectoryTrees source checkout also contains a package-native validation
runner and the complete packaged corpora under `tools/validation/` and
`inst/extdata/`:

| Package-native validation population | Inventory |
|---|---:|
| Adversarial/generated corpus | 626 trees |
| Registered known-true artifacts | 129 artifacts: 20 PST, 39 PST+SS, 21 PST+TT+SS, and 49 original-metric references |
| Empirical Anolis | 500 maps |
| Deterministic simulated Anolis | 1,000 maps |

Run all 129 package-native known-true comparisons from the TrajectoryTrees
root:

```bash
Rscript tools/validation/run_validation.R \
  --scopes=known-true \
  --output=tools/validation/output/known_true
```

### D. Additional cross-version regression checks

V34/V35 equality is additional evidence that the rewrite preserved existing
results. V34 is not the independent correctness authority. These comparisons
are separate from the logic tests and by-hand-confirmed references above.

Run construction, input-immutability, and V34/V35 comparison across every large
population with:

```bash
Rscript tools/validation/run_validation.R \
  --scopes=corpus,anolis-empirical,anolis-simulated \
  --compare-v34=true \
  --output=tools/validation/output/full_parity
```

The preserved regression evidence records exact PST/SS/TT agreement for 471
ultrametric corpus trees, 500 empirical Anolis maps, and 1,000 deterministic
simulated Anolis maps. The package-native runner does **not** itself execute the
seven/eight-family assertions; those are executed by the V35 release barrier
above. `tools/validation/` is kept in the source checkout for release
validation and is excluded from the built package by `.Rbuildignore`.

## 2. Automatic public-interface tests

Run the complete automatic suite with:

```bash
Rscript -e 'devtools::test(".", reporter = "summary")'
```

These tests live under `tests/testthat/`. They are intentionally fast and check
the installed public API. They are secondary to the numerical correctness
contracts above.

| Test or support file | What it verifies |
|---|---|
| `helper-fixtures.R` | Constructs the real `simmap` fixture used by the public constructor, batch, and distribution tests. The fixture contains a phylogeny, named state-duration maps, tip states, and mapped-edge totals. |
| `test-public-constructors.R` | Calls `make_trajectory_objects()` and checks the public `trajectory_tree` interface: fields, attributes, input-shape preservation, naming, and class retention. It does not prove the numerical contents of those fields. |
| `test-public-batches.R` | Writes constructed trajectory objects to batches and reloads them, checking the manifest, object names, ordering, and preserved `trajectory_tree` class. |
| `test-public-distributions.R` | Reduces multiple verified trajectory objects into a stable `trajectory_distribution` with the documented metadata and P/T/SS/TT surfaces. |
| `test-validation-data.R` | Verifies filtering and loading of the package validation-tree inventory by case label and tip count. |

## 3. Similarity-matrix and summary validation

Similarity validation begins only after the `trajectory_obj` construction path
above has been exercised. The 33-history and 500-tree comparisons construct one
`trajectory_obj` per history and pass that exact same object independently to
V2.1 and `TrajectoryTrees::trajectory_similarity()`.

| Test or support file | What it verifies |
|---|---|
| `test-public-similarity.R` | Checks the exported `trajectory_similarity()` interface, registered package-native C++ symbol, output structure, and argument validation. |
| `helper-scenario-mats-v21.R` | Loads the test-only V2.1 numerical authority and the three-tip reference constructors. |
| `reference/ScenarioMatsV2_1.R` | Supplies the readable plain-R numerical authority used only by tests. |
| `reference/construct_three_tip_reference_trees.R` | Constructs the 33 mapped-history reference cases used for exact differential comparison. |
| `test-similarity-matrix-summary-contract.R` | First checks hand-derived pairwise matrix entries. It then independently reconstructs synchronous and asynchronous tree-wide totals, proportions, means, and available-similarity totals from the returned matrices for ultrametric and non-ultrametric trees under both asynchronous metrics. |
| `test-similarity-v21-reference.R` | Compares the package-native matrices and summaries against V2.1 on all 33 hand-constructed histories, depths 0-3, both asynchronous metrics, and both return modes. It also checks the separately aggregated similarity vectors. |
| `test-similarity-validation-bundle.R` | Ensures the preserved 500-tree audit manifest and all 500 package-local RDS fixtures remain present, unique, and self-contained. It does not run the expensive audit. |

## Manual 500-tree V2.1/package-native audit

The preserved 500-tree audit lives at:

```text
tests/validation/scenariomats_v3_v21_500/
```

It is intentionally excluded from the ordinary package test command. Run the
complete audit, with no trailing arguments, from the package root:

```bash
Rscript tests/validation/scenariomats_v3_v21_500/run_audit.R
```

The runner applies three distinct validation layers to all 500 histories:

1. compare every package-native matrix against the V2.1 numerical authority;
2. check the independent matrix invariants and pairwise closure rules;
3. reconstruct tree-wide totals, proportions, means, and non-ultrametric
   available-similarity totals directly from those matrices and compare them
   with the reported summaries.

It uses transition depths 0-3, the Bhattacharyya and minimum asynchronous
metrics, and the complete and matrices-and-summaries return modes. It writes a
new timestamped directory under
`tests/validation/scenariomats_v3_v21_500/results/` and never overwrites a
previous run.

See `tests/validation/scenariomats_v3_v21_500/README.md` for the corpus and
output contracts.
