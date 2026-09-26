# owea 0.5.0

* c-optimality (one linear combination c'theta: a one-row `wb` or `grad_g`, or
  a one-parameter `subset`) is now certified by Elfving's bound instead of the
  sensitivity function. For any h with h'c = 1, c'M(w)^-c >= 1/(h'infor0 h +
  max_i h'A_i h) for every design w (Lagrangian dual of the convex criterion;
  strong duality holds), so the best h gives the c-optimal value and any h a
  valid lower bound. A design is certified when its value meets the bound
  within `eps0`, and `efficiency_lower_bound` is bound / value. This also
  certifies SINGULAR c-optimal designs (fewer support points than parameters),
  which the sensitivity function -- built from a pseudo-inverse of the
  information matrix -- reported as "did NOT converge" (e.g. estimating
  theta0 - theta1 = eta(-1) in a two-parameter logistic model: the one-point
  design at x = -1 is optimal, with value 4). The bound is exact on a grid or
  candidate set via a small linear program when `lpSolve` is installed (new in
  Suggests), otherwise (matrix-mode information, an existing design, no
  lpSolve) a derivative-free minimisation that can only be looser, never
  wrong; on the continuous path the maximum in the bound is over the points
  examined. New result fields `elfving_bound`, `elfving_gap`, `elfving_h` in
  `optimal_design()` and `verify_optimality()`; `exact_design()` inherits the
  reference design's certified bound.

* New `continuous = TRUE` option of `optimal_design()` (and `exact_design()`):
  the continuous covariates are searched in the CONTINUOUS design region
  instead of on a grid, so no `step_sequence` is needed and the cost no longer
  grows exponentially with the number of continuous covariates. The solver
  starts from the OWEA solution on a small coarse grid (`init_levels` = 3
  levels per continuous covariate times all factor levels, or user-supplied
  `init_points`), then alternates: polishing of the support-point locations
  (with the weights fixed, the criterion is minimised over all continuous
  coordinates by L-BFGS-B, using the closed-form gradient
  `-(w_i / v) d tr(P I(x)) / dx` at `x_i`, `P` being the matrix behind the
  directional derivative); the package's own Newton weight step (active-set
  engine) on the current support; merging of points closer than `merge_tol`;
  and adding the maximiser of the directional derivative over the region,
  found by multi-start L-BFGS-B (the support points, the box vertices and
  `n_starts` random starts, for every level combination of the factors). It
  stops when that maximum is at most `eps0` and a random audit of `n_audit`
  points (default 20,000; its best points polished locally) finds no
  violation. Factor covariates keep their levels. Works for D- and
  A-optimality, any `wb` / `subset` / `grad_g`, and an existing design. The
  certificate (`max_d`, `efficiency_lower_bound`) is the largest directional
  derivative over the points examined, not an exhaustive scan; `check_global`
  still runs a grid scan on top when `global_step` is given. The result
  carries `method = "continuous"`, `iterations`, `history`, `maximiser`,
  `n_audit`, `audit_max_d` and `jacobian`.
* Derivatives of the per-point information with respect to the continuous
  covariates: formula-style models (`link` + terms) now carry an analytic
  Jacobian and a vectorised evaluator (attributes `"jacobian"` and
  `"vectorized"` of `model_info_vector()`'s function); a user-supplied
  `info_vector` / `info_matrix` is differentiated by finite differences, or
  you can pass your own `info_jacobian`.
* `verify_optimality(continuous = TRUE)`: the grid-free counterpart of the
  check -- the maximum sensitivity over the continuous region by multi-start
  search plus the random audit (no `step` needed); returns `method`,
  `n_audit` and `audit_max_d`.
* App: when a covariate is continuous, the model step offers "How to search
  the continuous covariates": on a grid (enter the grid steps, as before) or
  directly in the continuous region (no steps). The continuous search runs
  without a random audit; its own check is the multi-start search of the
  region, and the verify panel offers the audit afterwards.
* App: the verify panel of a single-criterion approximate design checks the
  equivalence theorem either on a grid of the design box at step(s) you type
  (prefilled with the finest step of a grid computation, or one fortieth of
  each range after a continuous search; the usual large-grid safeguard) or by
  a random audit of the design region with the number of points you choose
  (default 20,000; the best points are polished and a multi-start search from
  the support points and box corners is added). After a grid computation the
  grid option notes that the computation's own grid was already verified, so a
  different step avoids repeating it. The result names the space it was
  checked over.
* App: the results page ends with an "R code for this analysis" section and a
  download button: a runnable script with the exact `optimal_design()` /
  `exact_design()` call the app made and, for approximate designs, the
  `verify_optimality()` check at the grid step entered in the verify panel;
  for an exact design whose simulation study has been run, the script also
  repeats the study with `simulate_design()` -- the exact design and every
  design it was compared with (the simple random sample, a custom design),
  pooled with any first stage -- and tabulates the mean squared errors. The
  compound branch has the same section: the objectives as a `components`
  list, the `compound_design()` / `compound_exact_design()` call, the
  `compound_criterion()` scoring of an approximate design (reference values
  reused), and the per-objective simulation study of an exact design.
* `simulate_design()` gains `existing` and `obs`: a first stage pooled with the
  design, as the app's simulation study does for a multistage design -- an
  existing design whose responses are simulated, or an observed data set whose
  responses are kept.
* Bug fix (app): after a continuous search, the simulation study's simple
  random sample failed with "each grid step must be a scalar ..." because it
  drew the runs from a grid the continuous search does not have. The sample is
  now drawn uniformly from the continuous region (factors uniformly over their
  levels); after a grid computation it is drawn from the finest-step grid.
* App: the simulation study's simple random sample from the finest-step grid
  no longer builds that grid, which grows exponentially with the number of
  covariates as the step shrinks (it used to be materialised and could exceed
  memory). The grid is a Cartesian product, so a grid index is drawn
  independently for each covariate: exactly uniform over the grid points,
  instant for any step. Both the single-criterion and the compound branch use
  it.
* Benchmarks (logistic model, D- and A-optimality): three covariates on
  [-2,2] x [-1,1] x [-3,3] -- the continuous search reaches the grid path's
  designs (marginally better criteria) in 0.3-0.5 s versus 2.4 s for the
  grid path with steps 0.1, 0.01 and a 0.05 global check; six covariates on
  [-2,2]^6, where a 0.05 grid would have 81^6 points -- 1 to 3.5 s.
* The default `continuous = FALSE` path is unchanged.

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
* Bug fix: `compound_exact_design()` accepts `reference_bound` (passed on to
  `compound_design()`). The app reuses the reference values of a previous run
  together with their bounds, and an exact compound design on such a rerun
  failed with "unused argument (reference_bound = ...)". The app's compound
  verify step now passes the bounds through as well.
* Clearer wording for the compound criterion value. The app's status line said
  "weighted average efficiency = 0.945 out of 1", which read as if 1 were
  attainable. It now explains that the value is the weighted average of the
  objectives' efficiencies, that it is the highest any single design can reach
  for these objectives and weights (certified by the max sensitivity), and that
  it would equal 1 only if one design were optimal for every objective at once.
  The print methods of `compound_design()` and `compound_exact_design()` say the
  same.

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
