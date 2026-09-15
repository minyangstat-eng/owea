# owea — Optimal Weights Exchange Algorithm for approximate optimal designs

`owea` finds **D- and A-optimal** experimental designs for linear, nonlinear, and
generalized linear models using the Optimal Weights Exchange Algorithm (OWEA) of
Yang, Biedermann & Tang (2013, *JASA* 108(504), 1411–1420). An alternative
**multiplicative-algorithm solver** (`solver = "MA"`, Yu 2010) is also available — see §9.

> ### ▶ Try the web app — no installation
> **<https://owea.shinyapps.io/owea-designs/>**
> Open it in any browser (works on phones and tablets too): describe your model
> with menus, click *Compute*, and get the design. No R, no setup.

The main entry point is **`optimal_design()`**. The rest of this README focuses
on it.

---

## What's new in 0.4.1

- **Efficiencies are guaranteed lower bounds** — every efficiency the package
  reports is now a single number, `efficiency_lower_bound`, certified relative to
  the *true* optimum (Theorems 4.5 and 4.6 of Becker & Yang), so a reference
  design that stopped short of the optimum can no longer overstate it. See the
  "Efficiency lower bounds" bullet in §8. Printed with four decimals.
- **App** — the model steps remind users that models or quantities of interest
  beyond the built-in ones can be used in the R package directly.

## What's new in 0.4.0

- **Active-set engine** — `optimal_design(engine = "active-set")` replaces the
  weight step of the OWEA solver by an active-set Newton method: the step is cut
  at the first weight that would reach zero and accepted under an Armijo decrease
  of the criterion, every point the step pushes to zero is pruned in one batch
  (the classic engine removes one point per Newton solve), and `add_per_iter`
  violating candidates are added per exchange iteration. It works for every
  criterion and quantity of interest (`wb`, `subset`, `grad_g`) and reaches the
  same designs as the classic engine. See §8.
- **Benchmarks** (`owea_engine_benchmark.R` in the repository root): on a
  401,841-point logistic grid with 4 parameters (D- and A-optimality, all
  parameters and a subset) both engines reach the same designs and solve in
  0.3–0.4 s. The new engine pays off when the weight step dominates — a large
  starting support or a design with hundreds of support points — where it cut
  solve times two- to four-fold in our tests. The default `engine = "classic"`
  is unchanged and the test suite passes identically.

## 0. Web app (no coding required)

For practitioners who would rather not write R, `owea` ships a point-and-click
**web app**. You describe the model with plain inputs — a model family, named
covariates (continuous ranges or factor levels), interactions picked by name, a
criterion, and assumed parameter values — and get the design as a table, a plot,
and the information matrix, with a CSV download. It hides the package's internal
conventions (factor coding, term indices, interaction codes).

The first screen asks which kind of problem you have:

- **One criterion (classical).** One model, one objective — a D- or A-optimal
  design, optionally for a subset of the parameters. This is the wizard
  described below, unchanged.
- **Several criteria at once (compound).** More than one objective, and they
  may come from **different models**. One design is found that serves them
  all, weighted as you choose — approximate or exact, with its own efficiency
  report, design scoring and simulation study. See §16.

Highlights of the wizard:

- **Grid step(s) / step sequence.** Each continuous covariate takes one grid
  step (`0.1`) or a comma-separated **step sequence** from coarse to fine
  (`0.5, 0.1, 0.02`): the whole range is searched at the coarsest step, then
  refined locally at each finer one — fine resolution without a huge grid.
- **Large-candidate-set safeguard.** The number of candidate design points is
  checked as soon as the model step is complete (and reported again on the
  Review step). Past 1,000,000 points you choose: adjust the grid step, switch
  to a step sequence, or proceed anyway.
- **Verify optimality of a design.** Paste or upload any design (e.g. an edited
  CSV) and check it under the *original* criterion over the design box at the
  finest grid step. If that grid would exceed 1,000,000 points you choose:
  cancel, run the full check anyway, or get the **criterion value only**.
- **Efficiency under a different criterion.** The computed design is evaluated
  against the design that is optimal for the other criterion, re-solved with
  the *exact same step sequence* as the original.
- **Existing designs and data.** An existing design or raw data set can be
  reused as a first stage; a data set can also supply the assumed parameter
  values by fitting the model (`fit_design()`).
- **Simulation study** (exact designs, in **both** branches). Data are
  generated at the design, the same model is fitted back, and the
  per-parameter **mean and median squared error** are averaged over the
  replications. The computed design can be compared against a simple random
  sample of the same size and against any design you paste in; every
  comparison design also gets a **ratio column** — its error divided by the
  computed design's, so above 1 the computed design is the more precise one
  for that parameter. The raw errors sit on each parameter's own scale, so the
  ratio is what compares designs across a row. Exact designs offer `n = 200`
  runs by default.

- **Run it locally** (needs the `shiny` and `DT` packages):

  ```r
  install.packages(c("shiny", "DT"))
  owea::run_owea_app()
  ```

- **Use it hosted — nothing to install:**
  **<https://owea.shinyapps.io/owea-designs/>**
  Open it in any browser (phones and tablets included) — no R, no Rtools, no
  setup. This is the easiest way to try the tool or share it with collaborators.

### Hosting on shinyapps.io (free)

shinyapps.io rebuilds the app's packages in the cloud. It fetches CRAN packages
automatically, but `owea` is **not on CRAN** and has a compiled C++ core, so the
one extra step is to serve `owea` from a **public GitHub repo** and install it
from there (that records the source the server rebuilds from). You need two free
accounts — **GitHub** (where the package source lives) and **shinyapps.io**
(where the app runs). The full script is [`deploy/deploy.R`](../deploy/deploy.R).

**One-time setup:**

1. Put the package on a **public GitHub repo** — push the contents of the
   `owea/` folder (so `DESCRIPTION` is at the repo root) to e.g.
   `github.com/<you>/owea`.
2. Install it locally **from GitHub** so `rsconnect` records where to rebuild it:
   ```r
   install.packages(c("remotes", "rsconnect"))
   remotes::install_github("<you>/owea")
   ```
3. Create a free account at <https://www.shinyapps.io>, then **Account → Tokens**
   and run the `rsconnect::setAccountInfo(...)` snippet it shows you.

**Deploy** (repeat this to publish updates):
```r
rsconnect::deployApp(
  appDir  = system.file("shiny", "owea-app", package = "owea"),
  appName = "owea-designs")
```
The first deploy compiles `owea` on the server (a few minutes), then prints a
public URL like `https://<you>.shinyapps.io/owea-designs/`. Free tier: 5 apps and
limited active hours/month. An institutional **Posit Connect** server uses the
same `deployApp()` call.

> WebAssembly/shinylive is **not** an option — RcppArmadillo does not compile to
> WASM. If you want a single-account, self-contained alternative, a Docker image
> (e.g. a free Hugging Face Space) can bundle `owea` directly, no GitHub needed.

The app is a thin wrapper over `optimal_design()` / `exact_design()`; everything
below documents the underlying R API it drives.

---

## 1. Quick start

```r
library(owea)

# logistic GLM with 3 factors; information vector f(x, theta),
# the per-point information matrix is f f'.
info_vec <- function(x, theta) {
  q   <- c(1, x[1], x[2], x[3])
  eta <- sum(q * theta)
  (exp(eta / 2) / (1 + exp(eta))) * q          # length-k vector
}

res <- optimal_design(
  info_vector   = info_vec,
  theta         = c(1, -0.5, 0.5, 1),
  design_box    = list(c(-2, 2), c(-1, 1), c(-3, 3)),
  step_sequence = c(0.2, 0.1, 0.05),
  p             = 0)                            # 0 = D, 1 = A, 2,3,... -> E

print_result(res)
res$support     # optimal support points (one per row)
res$weights     # their weights (sum to 1)
res$criterion   # criterion value (smaller is better)
res$max_d       # max directional derivative; <= eps0 confirms optimality
```

---

## 2. Specifying the model

Supply **exactly one** of `info_vector` or `info_matrix`. Each may be written
as `function(x)` or `function(x, theta)`:

| Argument | Returns | Per-point information |
|----------|---------|------------------------|
| `info_vector` | length-`k` vector `f` | `I_x = f f'` (rank one; the fast path) |
| `info_matrix` | `k × k` matrix | `I_x` directly (any model, incl. rank > 1) |

- Use `info_vector` when the per-point information is `f f'` — GLMs, nonlinear
  models with normal errors. It is faster.
- Use `info_matrix` when the per-point information is a general `k × k` matrix.

`x` is the covariate vector of one design point; index it as `x[1]`, `x[2]`, ….
If you write the function as `function(x, theta)`, supply `theta`; if you write
it as `function(x)` (with the parameters captured inside), `theta` is optional.

```r
# information matrix, two-argument form
info_mat <- function(x, theta) {
  q   <- c(1, x[1], x[2], x[3]); eta <- sum(q * theta)
  (exp(eta) / (1 + exp(eta))^2) * tcrossprod(q)     # nu(eta) * q q'
}
optimal_design(info_matrix = info_mat, theta = c(1, -0.5, 0.5, 1),
               design_box = list(c(-2,2), c(-1,1), c(-3,3)),
               step_sequence = c(0.2, 0.1, 0.05), p = 0)
```

---

## 3. Where the candidate points come from

Give **either** a continuous design region **or** a fixed candidate set:

- **Continuous region** — `design_box` (a list of `c(lo, hi)` pairs, one per
  covariate) and `step_sequence` (grid steps, coarsest first). The algorithm
  solves on the coarse grid, then refines only small neighbourhoods of the
  current support at each finer step, so the cost is dominated by the first
  stage. This reaches a fine resolution without ever materializing a fine grid.

  ```r
  optimal_design(info_vector = info_vec, theta = th,
                 design_box = list(c(-2,2), c(-1,1), c(-3,3)),
                 step_sequence = c(0.2, 0.1, 0.05, 0.02), p = 0)
  ```

- **Fixed candidate set** — `candidate_set`, an `n × N` matrix with one row per
  candidate point (it need not be a regular grid). When given,
  `design_box`/`step_sequence` are ignored.

  ```r
  X <- candidate_grid(list(c(0, 3)), step = 0.05)   # or make_grid(), or your own matrix
  optimal_design(info_vector = info_vec, theta = th, candidate_set = X, p = 0)
  ```

`candidate_grid(design_box, step)` and `make_grid(lower, upper, by)` build a
rectangular grid; or pass any `n × N` numeric matrix of your own points.

---

## 4. Optimality criterion `p`

Only **D-** and **A-optimality** are supported; any other `p` raises an error.

| `p` | Criterion | Reported value |
|-----|-----------|----------------|
| `0` | D-optimal | `log|Σ| / v` |
| `1` | A-optimal | `tr(Σ) / v` |

where `Σ = (∂g/∂θ) M⁻¹ (∂g/∂θ)'`, `M` is the design information matrix, and `v`
is the number of quantities of interest. Smaller `criterion` is better. `p` must
be `0` (D) or `1` (A).

---

## 5. Quantity of interest (full vector, a subset, or a function of θ)

By default the design targets the **full** parameter vector `θ`. To target
something else, pass **at most one** of:

| Argument | For | Form |
|----------|-----|------|
| `subset` | a **subset** of the parameters | integer indices, e.g. `subset = c(2, 4)` |
| `grad_g` | a general differentiable `g(θ)` | a function `theta -> v × k` Jacobian `∂g/∂θ` |
| `wb` | a general **linear** `g(θ) = W θ` | a constant `v × k` matrix `∂g/∂θ` |

`subset = c(2, 4)` is the convenient form of a row-selecting Jacobian; `grad_g`
and `wb` are for any other (possibly non-subset) quantity. For example, interest
in `g(θ) = (θ₁, θ₂ − θ₄)`:

```r
grad_g <- function(theta) matrix(c(1,0,0,0,
                                   0,1,0,-1), nrow = 2, byrow = TRUE)
optimal_design(info_matrix = info_mat, theta = th, design_box = box,
               step_sequence = steps, p = 1, grad_g = grad_g)
```

> Note: for a partial-parameter / general `g(θ)` the optimum can drive the
> *full* information matrix toward singularity, where the equivalence theorem
> cannot certify optimality. The run does not crash (a pseudo-inverse is used),
> but may end with `converged = FALSE` while the criterion has stabilized.

---

## 6. Augmenting an existing design (multistage)

To allocate new runs so the **combined** design (an existing design ξ₀ plus the
new one) is optimal, pass `xi0_points`, `xi0_weights`, and the sample sizes
`n0` (existing) and `n1` (new):

```r
optimal_design(info_vector = info_vec, theta = th, design_box = box,
               step_sequence = steps, p = 0,
               xi0_points  = matrix(c(0, 1, 2, 3), ncol = 1),
               xi0_weights = rep(0.25, 4), n0 = 40, n1 = 80)
```

`n0 = 0` (the default) is a single-stage / locally optimal design.

---

## 7. Optional merging of neighbouring support points

On a discrete grid the optimum can be represented as two adjacent points
splitting one support point's weight. Set `merge = TRUE` to merge support points
that are close together (weighted centroid, then re-optimise the weights):

- For a `design_box`, the tolerance is `merge_factor * step` at each stage
  (`merge_factor` default `1.5`).
- For a `candidate_set`, it defaults to `merge_factor` times the smallest
  positive per-coordinate gap; override with `merge_atol`.

`merge = FALSE` (the default) returns the design exactly as the engine found it.

---

## 8. Other arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `theta` | `NULL` | parameter values; required only for a `function(x, theta)` model, `grad_g`, or an existing design |
| `solver` | `"owea"` | design algorithm: `"owea"` (the exchange engine) or `"MA"` (the multiplicative algorithm as a direct solver, see §9) |
| `init_method` | `"auto"` | starting support for the OWEA engine: `"auto"` (IBOSS for vector input, minmax for matrix input), `"minmax"`, `"minmaxmedian"`, `"random"`, `"iboss"`, or `"MA"` (a multiplicative-algorithm warm start — keeps the `k+1` highest-weight points; often fewer engine iterations) |
| `ma_max_iter` | `100` | multiplicative-algorithm iteration cap. For `init_method = "MA"` it applies exactly; for `solver = "MA"` the effective cap is `max(ma_max_iter, 10000)` |
| `auto_warm_start` | `TRUE` | if a `candidate_set` solve fails from a cold start, automatically retry warm-started from a quick coarse multistage solve (OWEA engine only) |
| `check_global` | `FALSE` | (`design_box` path) after converging, verify the design over a fine grid spanning the whole box; reports `global_max_d` / `global_check` |
| `global_step` | `NULL` | grid step for the `check_global` verification (default the finest `step_sequence` step) |
| `max_iter` | `100` | maximum outer iterations per stage |
| `eps0` | `1e-6` | stopping threshold on the directional derivative |
| `engine` | `"classic"` | weight step of the OWEA solver: `"classic"` (the original damped Newton, one point pruned per Newton solve) or `"active-set"` (ratio-test/Armijo Newton step, batch pruning, batched additions — see below) |
| `add_per_iter` | `1` | `engine = "active-set"` only: number of violating candidates added per exchange iteration (at most 20% of the current support) |
| `accept_tol` | `1e-9` | a refinement stage is kept only if it converges and does not worsen the criterion by more than this |
| `verbose` | `FALSE` | print stage-by-stage progress |

### The active-set engine (`engine = "active-set"`)

Both engines run the same exchange loop — add the candidate with the largest
directional derivative, re-optimise the weights on the support, repeat until
`max_d <= eps0` — and converge to the same design. They differ in the weight
step:

- **classic** — a damped Newton step whose damping is halved whenever a trial
  leaves the simplex; after each converged solve the single smallest-weight
  point is removed and the solve repeated. Pruning a large starting support
  therefore costs one Newton solve per removed point.
- **active-set** — the Newton direction is cut at the first weight that would
  reach zero and accepted only under an Armijo decrease of the criterion; all
  points the step pushes to zero are removed in one batch. A point is never
  dropped while its directional derivative is positive (it is under-weighted,
  not superfluous), the support never falls below the rank-aware minimum, and a
  batch is rejected — falling back to a single removal — if it would make the
  quantity of interest non-estimable. With `add_per_iter > 1` the worst
  violators are added several at a time.

Use it when the starting support is large (the default minmax start hands
`2 × ncol(candidate_set)` points to the weight step), when the optimal design
has hundreds of support points, or when a solve seems to spend its time in the
weight step. For large candidate sets a good combination is

```r
optimal_design(info_matrix = f, candidate_set = X, p = 0,
               engine = "active-set", add_per_iter = 5,
               init_method = "MA", ma_max_iter = 20, max_iter = 1000)
```

Raise `max_iter` for problems with many support points: each exchange
iteration adds at most `add_per_iter` points (one for the classic engine), so a
design with 300 support points needs at least that many iterations from a small
start.

### Return value

A list with `support` (an `n × N` matrix, one design point per row), `weights`,
`criterion`, `max_d` (max directional derivative; `<= eps0` confirms optimality),
`converged`, `times`, `grid_sizes`, `total_time`, `box_lo`, `box_hi`, `p`, and
`global_max_d` / `global_check` (`NA` unless `check_global = TRUE` or a
`candidate_set` was used).

### Convergence and global optimality

Always check `res$converged` / `res$max_d`. The package also warns you:

- **Not converged** — `optimal_design()` / `owea()` warn when the returned design
  is not optimal (`max_d > eps0`); raise `max_iter`/`max_outer`, coarsen the grid,
  or warm-start.
- **Local (multistage) convergence** — on the `design_box` path a converged
  design is optimal only over the *refined neighbourhood grids*, not the whole
  box; a warning says so. Pass `check_global = TRUE` to verify over a fine grid
  spanning the whole box (`global_check` is then `TRUE`/`FALSE`).
- **Hard problems** — if a cold `candidate_set` solve stalls (common for many
  parameters on a large grid), `auto_warm_start = TRUE` (default) retries it
  warm-started from a coarse multistage solve, usually reaching the global
  optimum transparently. That coarse grid spans the bounding box of the
  candidate columns, so it is only feasible for candidate sets with a handful
  of columns; for candidate sets with many columns use `init_method = "MA"`
  with a larger `max_iter` and `engine = "active-set"` instead.
- **Efficiency lower bounds** — a design that stops at `max_d > eps0` is not
  exactly optimal, and an efficiency measured against such a design could be
  overstated. So the package reports one number, `efficiency_lower_bound`: a
  guaranteed lower bound on the efficiency relative to the *true* optimum over
  the design space checked. For a single design it is computed from the design
  alone (Theorems 4.5 and 4.6 of Becker & Yang, *Post Hoc Control Group
  Selection via Constrained Optimal Design*): `1 / (1 + max_d / v)` for
  D-optimality and `1 - max_d / criterion` for A-optimality, with `v` the number
  of parameters of interest; it equals 1 at a certified optimum.
  `optimal_design()` and `verify_optimality()` return it directly. Where a
  design is compared with a derived reference design — `exact_design()`,
  `compound_design()`'s per-component values and cross-efficiency table, the
  app's efficiency panels — the ratio is multiplied by the reference's own
  bound, so it stays valid even when the reference solve stopped short.
  Bounds are printed with four decimals (e.g. `efficiency >= 99.9996%`).

---

## 9. An alternative solver — the multiplicative algorithm (`solver = "MA"`)

By default `optimal_design()` uses the OWEA exchange engine (`solver = "owea"`).
Passing `solver = "MA"` (alias `"multiplicative"`) instead solves the approximate
design with the **multiplicative algorithm** (Yu 2010) run to convergence, and
returns *its* design directly — no exchange engine on top.

```r
optimal_design(info_matrix = info_mat, theta = th,
               design_box = list(c(-1,1), c(-1,1), c(-1,1), c(-1,1), c(-1,1)),
               step_sequence = c(1), p = 0, solver = "MA")
```

- **Coverage.** It handles D- and A-optimality (`p = 0`/`1`), any quantity of
  interest (`subset` / `grad_g` / `wb`), and existing designs (`n0 > 0`). For a
  criterion outside its guaranteed-convergent class (c-optimality / rank-1 `wb`),
  each per-grid solve that cannot certify optimality **falls back automatically**
  to the OWEA engine, so the result is always correct.
- **Speed.** `"owea"` is generally faster; but for a large dimension of the
  information matrix (many parameters / a big grid) `"MA"` has the advantage. It
  returns the multiplicative design, whose support may be larger (weight spread
  over more points) than the sparse exchange-engine one.
- **`ma_max_iter`** (default `100`) caps the iterations; as the solver the
  effective cap is `max(ma_max_iter, 10000)` — raise it above `10000` for a
  slowly converging problem.

When `solver = "MA"`, `init_method` and `auto_warm_start` are unused.
`exact_design()` also accepts `solver = "MA"` (it solves the reference
approximate design with MA, then rounds and exchanges as usual).

> **`init_method = "MA"` is different.** That keeps the OWEA engine but *warm-starts*
> it from a quick multiplicative pass (often fewer engine iterations); `solver = "MA"`
> *replaces* the engine. They are independent knobs.

---

## 10. Exact designs for a given sample size — `exact_design()`

`optimal_design()` returns an **approximate** design (continuous weights summing
to 1). When you must run an experiment with a fixed number of runs `n`, use
**`exact_design()`** to get an **exact** design: an integer allocation
`n₁, …, n_m` of the `n` runs over design points (`Σ nᵢ = n`).

`exact_design()` takes the same model/region arguments as `optimal_design()`
plus the sample size `n` (first argument). It computes the approximate optimum
internally, then (i) rounds the weights to integer counts by largest-remainder
apportionment, (ii) repairs the total to exactly `n` using the sensitivity
(directional-derivative) function — peeling a run off the most *inefficient*
support point or adding one to the most *informative* candidate — and (iii)
improves the design by random exchanges. It also reports the design's
**efficiency** relative to the approximate optimum (in `(0, 1]`, `1` = optimal),
with the exact design evaluated as an approximate design with weights `nᵢ / n`.

```r
r <- exact_design(
  n             = 30,
  info_vector   = info_vec, theta = c(1, -0.5, 0.5, 1),
  design_box    = list(c(-2, 2), c(-1, 1), c(-3, 3)),
  step_sequence = c(0.2, 0.1, 0.05),
  p             = 0, seed = 1)

print(r)
r$counts        # integer runs per support point (sum to n)
r$support       # the exact support points (one per row)
r$efficiency    # guaranteed LOWER BOUND on the exact design's efficiency
r$criterion     # the exact design's criterion value
```

The reported efficiency is a **lower bound**: it compares the exact design to the
approximate optimum, which is at least as good as any exact design, so the true
efficiency is *at least* this value — hence the printout shows `efficiency >= …`.

With a **step sequence**, the exchanges in (iii) search a *coarse* grid over the
whole box (the coarsest step in the sequence) together with a *fine*
neighbourhood of the approximate support at the finest step — global reach where
it is needed, fine resolution where it matters. The finest grid is never built
over the whole region: its size grows as the reciprocal of the final step raised
to the number of continuous covariates, so a fine last stage would otherwise
cost far more here than it does in the approximate solve. A single-step
sequence is unaffected (coarse and fine coincide, and the candidate set is the
full grid), and `candidate_set` overrides the construction entirely.

Useful arguments: `max_exchange` (number of random exchanges, default `1000`),
`seed` (reproducibility), `snap_support` (`design_box` path: snap the approximate
support to the finest grid, default `TRUE`; or append it off-grid),
`solver` / `ma_max_iter` (solve the reference approximate design with the
multiplicative algorithm — `solver = "MA"`, see §9), and
`check_global` / `global_step` / `global_max_points` (passed through to the
internal `optimal_design()` to certify the approximate reference over the whole
box — see §8). `print(res)` (or `print_result(res)`) reports the design's
efficiency lower bound as a percentage and its integer run counts. For a multistage call
(`xi0_*`, `n0 > 0`), `n` is the exact size of the **new** stage.
`exact_design()` errors if `n` is below the minimum support needed for a
non-singular information matrix.

---

## 11. Verifying a given design — `verify_optimality()`

To check whether **any** design (support points + weights) is optimal — not just
one the package computed — use `verify_optimality()`. It takes the same model
and quantity-of-interest arguments as `optimal_design()`, plus the design and
the design space (`design_box` + a single `step`, or a `candidate_set`), and
reports the **maximum sensitivity** over the space (0 at the optimum, by the
equivalence theorem), the criterion value, and the information matrix. The
design is validated first: weights must be nonnegative and sum to 1, and
support points off the grid are flagged with a warning (checked coordinate by
coordinate — instant even for a huge grid).

```r
v <- verify_optimality(res$support, res$weights, info_vector = info_vec,
                       theta = th, design_box = box, step = 0.05, p = 0)
v$is_optimal$value    # TRUE if max_sensitivity <= tol
v$max_sensitivity     # <= 0 (up to tol) at the optimum
v$criterion
```

The optimality check evaluates the sensitivity at **every** point of the
`design_box` + `step` grid. If that grid exceeds `max_points` (default `1e6`)
you are asked to **abort**, **proceed anyway**, or compute the
**criterion only** (`criterion_only = TRUE`, also available directly): the
design is still validated and its criterion and information matrix returned,
but no grid is built and optimality is not assessed (`max_sensitivity` is
`NA`).

---

## 12. Related functions

- `owea()` + `DesignProblem()` — a lower-level interface for a fixed candidate
  set (the same engine `optimal_design()` uses), returning the design directly
  from a `DesignProblem` object.
- `appro_opt()` / `appro_opt_seq()` — information-vector API taking a precomputed
  `k × N` matrix or a sequential search.
- `make_grid()` / `candidate_grid()` — build a rectangular candidate grid.
- `design_information()` / `infor_matrix()` — information matrix of a design.
- `find_best_point()` — equivalence-theorem check of a design on any grid.
- `fit_design()` — fit the model to a data set (used by the app to turn data
  into assumed parameter values and a first-stage design).

---

## 13. Installation

The package compiles C++ on install, so every user needs:

1. **R ≥ 4.0**.
2. A **C++ toolchain**:
   - **Windows** — **Rtools** matching your R version
     (<https://cran.r-project.org/bin/windows/Rtools/>; e.g. Rtools43 for R 4.3.x).
   - **macOS** — `xcode-select --install`.
   - **Linux** — `g++` / `build-essential`.
3. The R packages **Rcpp** and **RcppArmadillo** (installed automatically by the
   methods below, or `install.packages(c("Rcpp", "RcppArmadillo"))`).

Check the toolchain: `pkgbuild::has_build_tools(debug = TRUE)` should be `TRUE`.

```r
# from GitHub:
remotes::install_github("minyangstat-eng/owea")
# or from a source folder / tarball:
install.packages("/path/to/owea", repos = NULL, type = "source")
```

Then `library(owea)`.

---

## 14. Demos and examples

```r
demo(package = "owea")                # list demos
demo("demo_unified")                  # vector vs matrix input, two-stage, optional merge
demo("demo_general")                  # multistage A-optimal logistic GLM
demo("demo_subset")                   # subset / grad_g
demo("demo_exact")                    # exact (integer) designs for a given sample size

source(system.file("examples", "examples.R",           package = "owea"))
source(system.file("examples", "example_logistic3d.R", package = "owea"))
```

---

## 15. References

Yang, M., Biedermann, S. & Tang, E. (2013). On Optimal Designs for Nonlinear
Models: A General and Efficient Algorithm. *JASA* 108(504), 1411–1420.

Yu, Y. (2010). Monotonic convergence of a general algorithm for computing optimal
designs. *Annals of Statistics* 38(3), 1593–1606. (The multiplicative algorithm
behind `solver = "MA"`.)

Becker, E. & Yang, M. Post Hoc Control Group Selection via Constrained Optimal
Design. Manuscript. (Theorems 4.5 and 4.6: the computable bounds on the
D- and A-optimality gaps behind every reported `efficiency_lower_bound`.)

---

## 16. Compound designs: one design, several objectives — `compound_design()`

Everything above finds a design for **one** criterion. When you have several
objectives at once — and especially when they come from **different models** —
use `compound_design()`. It maximises the weighted average of the component
**efficiencies**

```
Psi_alpha(xi) = sum_j  alpha_j * Psi_j(xi) / Psi_j*
```

where component *j* has its own model, its own quantity of interest, and its own
type, and `Psi_j*` is the value the design optimal for component *j* **alone**
would achieve. Each component criterion is stated *per parameter*,

| type | `p` | `Psi_j` |
|------|-----|---------|
| D | `0` | `det(S_j^-1)^(1/v_j)` |
| A | `1` | `v_j / tr(S_j)` |

with `S_j = G_j M_j^-1 G_j'`, so components of different dimension are
comparable and `Psi_alpha` lies in `(0, 1]` and reads directly as an average
efficiency. The design region is shared — one experiment — but the information
each model takes from it is not.

```r
f1 <- function(x, th) { q <- c(1, x[1], x[2]); e <- sum(q * th)
                        (exp(e / 2) / (1 + exp(e))) * q }          # logistic
f2 <- function(x, th) { q <- c(1, x[1], x[2], x[1] * x[2]); e <- sum(q * th)
                        (exp(e / 2) / (1 + exp(e))) * q }          # + interaction

res <- compound_design(
  components = list(
    list(info_vector = f1, theta = c(0.5, 1, -1),      p = 0, name = "main"),
    list(info_vector = f2, theta = c(0.5, 1, -1, 0.5), p = 0, name = "interaction")),
  alpha         = c(0.5, 0.5),
  design_box    = list(c(-2, 2), c(-2, 2)),
  step_sequence = c(0.5, 0.1))

print(res)
res$efficiency         # what each objective gets from this one design
res$cross_efficiency   # ... and what a single-objective design would have cost
```

`print()` reports the **cross-efficiency table**, which is the point of the
whole exercise: every design scored under every objective.

```
                              logistic main logistic + inter Poisson A(slopes)
optimal for logistic main             1.000            0.305             0.193
optimal for logistic + inter          0.697            1.000             0.760
optimal for Poisson A(slopes)         0.716            0.403             1.000
THIS DESIGN                           0.859            0.889             0.941
```

The D-optimal design for the main-effects model keeps only 19 % efficiency under
the Poisson component; the compound design is above 85 % under all three.

Components take the same model arguments as `optimal_design()` — either
`info_vector`/`info_matrix` + `theta`, or a formula-style spec (`link`, `f`,
`x`, `fx`, `xx`, `ff`, `ncat`, `coding`) — plus `p` (0 = D, 1 = A) and
optionally `subset` / `grad_g` / `wb`. Existing designs (`xi0_points`, `n0`,
`n1`) work as elsewhere; the runs are shared, the information is per model.

### Exact compound designs — `compound_exact_design()`

The compound counterpart of `exact_design()`: an integer allocation of `n` runs
that serves every objective. It takes the same arguments plus `n`, computes the
approximate compound optimum, apportions it to whole runs, and improves it by
random exchanges — the acceptance test reversed, since `Psi_alpha` is
*maximised*.

```r
ex <- compound_exact_design(
  n = 60,
  components    = list(
    list(info_vector = f1, theta = c(0.5, 1, -1),      p = 0, name = "main"),
    list(info_vector = f2, theta = c(0.5, 1, -1, 0.5), p = 0, name = "interaction")),
  design_box    = list(c(-2, 2), c(-2, 2)),
  step_sequence = c(0.5, 0.1),
  seed = 1)

ex$counts            # integer runs per support point (sum to n)
ex$efficiency_exact  # LOWER bound, vs the approximate compound optimum
ex$n_candidates      # how many candidate points the exchanges searched
```

The candidate set follows the same coarse-box-plus-fine-neighbourhood rule as
`exact_design()` (§10), so a fine last stage in the `step_sequence` costs little
here: the finest grid is never materialised over the whole region.

Related functions:

- `compound_criterion()` — score **any** design under the compound criterion
  and, given a design space, get the maximum sensitivity that certifies
  optimality. The scan visits every point of the `design_box` + `step` grid, so
  it is capped by `max_points` (default `1e6`) exactly as
  `verify_optimality()` is: past the cap you choose to abort, proceed anyway,
  or take `criterion_only = TRUE` — the design is still scored, but no grid is
  built and optimality is not assessed (`max_d` and `is_optimal` are absent).
- `compound_sensitivity()` — a thin wrapper when the sensitivity function over
  a candidate set is what you are after.
- Set `efficiency = FALSE` to weight the raw `Psi_j` instead. The two give
  different optima unless every `Psi_j*` is equal, and raw values of different
  models are not comparable, so the default is usually what you want.
- Supply `psi_star` to skip the reference solves. Computing `Psi_j*` costs one
  full solve **per component** before the compound solve even starts, which
  dominates a compound run; passing values from an earlier run with the same
  objectives, region and existing design reuses them exactly. (The web app
  does this for you: re-computing after changing only the weights `alpha`
  skips them, since `Psi_j*` does not depend on `alpha`.)

Fitting `theta` from a data set is the one classical-branch feature with no
compound counterpart; use `fit_design()` and pass the estimates in as each
component's `theta`.

The theory — the compound equivalence theorem, the weight gradient and
Hessian, and the derivation of both — is in `theory/compound-criterion.tex`
and `theory/second-derivative-derivation.tex`.
