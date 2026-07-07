# owea — Optimal Weights Exchange Algorithm for approximate optimal designs

`owea` finds **D- and A-optimal** experimental designs for linear, nonlinear, and
generalized linear models using the Optimal Weights Exchange Algorithm (OWEA) of
Yang, Biedermann & Tang (2013, *JASA* 108(504), 1411–1420).

> ### ▶ Try the web app — no installation
> **<https://owea.shinyapps.io/owea-designs/>**
> Open it in any browser (works on phones and tablets too): describe your model
> with menus, click *Compute*, and get the design. No R, no setup.

The main entry point is **`optimal_design()`**. The rest of this README focuses
on it.

---

## 0. Web app (no coding required)

For practitioners who would rather not write R, `owea` ships a point-and-click
**web app**. You describe the model with plain inputs — a model family, named
covariates (continuous ranges or factor levels), interactions picked by name, a
criterion, and assumed parameter values — and get the design as a table, a plot,
and the information matrix, with a CSV download. It hides the package's internal
conventions (factor coding, term indices, interaction codes).

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
| `init_method` | `"auto"` | starting support: `"auto"` (IBOSS for vector input, minmax for matrix input), `"minmax"`, `"minmaxmedian"`, `"random"`, `"iboss"` |
| `auto_warm_start` | `TRUE` | if a `candidate_set` solve fails from a cold start, automatically retry warm-started from a quick coarse multistage solve |
| `check_global` | `FALSE` | (`design_box` path) after converging, verify the design over a fine grid spanning the whole box; reports `global_max_d` / `global_check` |
| `global_step` | `NULL` | grid step for the `check_global` verification (default the finest `step_sequence` step) |
| `max_iter` | `100` | maximum outer iterations per stage |
| `eps0` | `1e-6` | stopping threshold on the directional derivative |
| `accept_tol` | `1e-9` | a refinement stage is kept only if it converges and does not worsen the criterion by more than this |
| `verbose` | `FALSE` | print stage-by-stage progress |

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
  optimum transparently.

---

## 9. Exact designs for a given sample size — `exact_design()`

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

Useful arguments: `max_exchange` (number of random exchanges, default `1000`),
`seed` (reproducibility), `snap_support` (`design_box` path: snap the approximate
support to the finest grid, default `TRUE`; or append it off-grid), and
`check_global` / `global_step` / `global_max_points` (passed through to the
internal `optimal_design()` to certify the approximate reference over the whole
box — see §8). `print(res)` (or `print_result(res)`) reports the design's
efficiency lower bound as a percentage and its integer run counts. For a multistage call
(`xi0_*`, `n0 > 0`), `n` is the exact size of the **new** stage.
`exact_design()` errors if `n` is below the minimum support needed for a
non-singular information matrix.

---

## 10. Related functions

- `owea()` + `DesignProblem()` — a lower-level interface for a fixed candidate
  set (the same engine `optimal_design()` uses), returning the design directly
  from a `DesignProblem` object.
- `appro_opt()` / `appro_opt_seq()` — information-vector API taking a precomputed
  `k × N` matrix or a sequential search.
- `make_grid()` / `candidate_grid()` — build a rectangular candidate grid.
- `design_information()` / `infor_matrix()` — information matrix of a design.
- `find_best_point()` — equivalence-theorem check of a design on any grid.

---

## 11. Installation

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
remotes::install_github("USERNAME/owea")
# or from a source folder / tarball:
install.packages("/path/to/owea", repos = NULL, type = "source")
```

Then `library(owea)`.

---

## 12. Demos and examples

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

## 13. Reference

Yang, M., Biedermann, S. & Tang, E. (2013). On Optimal Designs for Nonlinear
Models: A General and Efficient Algorithm. *JASA* 108(504), 1411–1420.
