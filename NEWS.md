# owea 0.4.1

* Efficiencies are now reported as guaranteed lower bounds relative to the TRUE
  optimum, in one field named `efficiency_lower_bound` (the former `efficiency`
  fields of `exact_design()`, `compound_design()`, `compound_criterion()` and
  `compound_exact_design()` are renamed; `efficiency_exact` becomes
  `efficiency_exact_lower_bound`). A design returned at `max_d > eps0` is not
  exactly optimal, so an efficiency measured against it could be overstated.
  `optimal_design()` and `verify_optimality()` gain `efficiency_lower_bound`,
  computed from the design alone: `1 / (1 + max_d / v)` for D and
  `1 - max_d / criterion` for A (Theorems 4.5 and 4.6 of Becker & Yang, *Post
  Hoc Control Group Selection via Constrained Optimal Design*; owea's `max_d` is
  their optimality gap). `exact_design()`'s value is the ratio to the
  approximate design times that design's bound (`approx$efficiency_bound`);
  `compound_design()`'s per-component values and `cross_efficiency` are the
  ratios to the reference designs times their bounds (`reference_bound`; also
  accepted as an argument alongside a supplied `psi_star`), and
  `criterion_bound` is the compound design's own bound. The app shows these
  values. Bounds are printed with four decimals (e.g. `efficiency >= 99.9996%`).
* App: the model step and the compound objectives step remind users that models
  or quantities of interest beyond the built-in ones can be used in the R
  package directly through `info_vector` / `info_matrix` and `wb` / `grad_g`.

# owea 0.4.0

* New `engine` argument of `optimal_design()`. `engine = "active-set"` replaces
  the weight step of the OWEA solver by an active-set Newton method: the step is
  cut at the first weight that would reach zero and accepted under an Armijo
  decrease of the criterion, every point the step pushes to zero is pruned in one
  batch (the classic engine removes one point per Newton solve), and
  `add_per_iter` violating candidates are added per exchange iteration. A point
  is never dropped while its directional derivative is positive, the support
  never falls below the rank-aware minimum, and a batch that would make the
  quantity of interest non-estimable falls back to a single removal. Applies to
  every criterion and to `wb` / `subset` / `grad_g`.
* New `add_per_iter` argument (active-set engine only; at most 20% of the current
  support is added per iteration).
* The default `engine = "classic"` is the previous code path, unchanged; the
  test suite passes identically. The C++ entry point `appro_opt_cpp()` gained the
  trailing arguments `mode` and `add_per_iter` with defaults that reproduce the
  old behaviour.
* Benchmarks: `owea_engine_benchmark.R` in the repository root. On a
  401,841-point logistic grid with 4 parameters (D- and A-optimality, all
  parameters and a subset) both engines reach the same designs and solve in
  0.3-0.4 s. The new engine pays off when the weight step dominates (a large
  starting support or a design with hundreds of support points), where it cut
  solve times two- to four-fold in our tests.
* README: new "What's new" section, the two arguments in the argument table, a
  section on the active-set engine, and a note that the automatic coarse-grid
  warm start is only feasible for candidate sets with few columns.
* Bug fix, `compound_design()`: a component's reference optimum `Psi*` could be
  far too small, which inflated every efficiency reported for that component
  above 1. The compound engine's single-component solve could prune the support
  into a singular information matrix (e.g. 6 points for a 7-parameter logistic
  model with interactions whose quantity of interest was the three main-effect
  slopes, A-optimality), after which the exchange loop cycled without
  converging. Two fixes: the engine never removes a point if any component's
  information matrix would become singular, and each reference value is now
  cross-checked against `optimal_design()` for that component, keeping the
  better design. `compound_design()` also warns explicitly when an efficiency
  exceeds 1.
* Reporting, `compound_design()` / `compound_criterion()` /
  `compound_exact_design()`: the per-component criterion values are now also
  returned on the scale `optimal_design()` reports (`component_criterion`,
  `component_criterion_star`: log det Sigma / v for D, tr(Sigma) / v for A,
  i.e. -log(Psi) and 1/Psi), and the print methods and the app's summary table
  show these instead of Psi, so the single-criterion and compound results agree
  number for number. The compound criterion itself, its optimisation and the
  efficiencies are unchanged; `psi` and `psi_star` remain in the results.
* App: when an exact design is chosen (single-criterion and compound flows), the
  design step now shows a note that the results page offers a simulation study
  for exact designs, so users know the design can be checked by simulating
  responses and re-estimating the parameters.
# owea 0.3.0

* Compound designs (`compound_design()`, `compound_exact_design()`), exact-design
  speedups and simulation MSE ratios.
* App: step sequences and candidate-set safeguards; `verify_optimality()` gained
  `criterion_only`.
