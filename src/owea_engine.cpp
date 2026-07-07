// =====================================================================
//  owea_engine.cpp -- unified C++ (RcppArmadillo) core for approximate
//                     optimal design (OWEA).
//
//  R port / merge of the two Julia-derived implementations of the
//  Optimal Weights Exchange Algorithm of Yang, Biedermann & Tang
//  (2013, JASA 108:1411-1420).
//
//  ONE engine handles BOTH per-point information representations:
//
//    info_mode = 0  (VECTOR / fast path)
//        info_data is k x N; column n is f(x_n, theta), so the per-point
//        information matrix is I_n = f f'.  The Newton step exploits this
//        rank-1 structure and the equivalence-theorem scan is the single
//        batched product  sum( (part * IVA) % IVA, 0 ).
//
//    info_mode = 1  (MATRIX / general path)
//        info_data is (k*k) x N; column n is vec(I_n), the full k x k
//        Fisher information matrix at x_n stored column-major.  This
//        covers models whose per-point information is not rank one.  The
//        scan becomes the single product  vec(part)' * I_mat.
//
//  EXISTING (MULTISTAGE) DESIGNS
//        The combined information of a two-stage design is
//            a * I_xi0  +  b * I_xi ,   a = n0/(n0+n1),  b = n1/(n0+n1).
//        The R layer passes  infor0 = a * I_xi0  and PRE-SCALES the
//        candidate information by b (vectors by sqrt(b), matrices by b),
//        so the assembled candidate information equals b * I_xi and
//        spd_inv(infor + infor0) = spd_inv(a I_xi0 + b I_xi).  With
//        infor0 = 0 (single stage) everything reduces to the original
//        single-stage code exactly.
//
//  Conventions
//   * Design-point indices are 1-based across the R boundary.
//   * pp = 0 -> D-optimality, pp = 1 -> A-optimality, pp >= 2 -> Phi_p.
//   * The reported criterion is the NORMALISED Phi_p value:
//        pp == 0 :  log|Sigma| / v
//        pp  > 0 : ( tr(Sigma^p) / v )^(1/p) ,  Sigma = wb M^{-1} wb'.
// =====================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------
//  spd_inv  --  robust inverse of a symmetric positive-(semi)definite
//               matrix (Cholesky first, then an increasing ridge, then
//               pinv as a last resort).  Kept large in singular
//               directions on purpose (see the original owea.cpp note).
// ---------------------------------------------------------------------
static mat spd_inv(const mat& A) {
    mat   S = 0.5 * (A + A.t());
    uword n = S.n_rows;
    mat   R;
    if (inv_sympd(R, S)) return R;
    mat    Id  = eye<mat>(n, n);
    double reg = 1e-12 * std::max(1.0, trace(S) / (double) n);
    for (int t = 0; t < 12; ++t) {
        if (inv_sympd(R, S + reg * Id)) return R;
        reg *= 10.0;
    }
    return pinv(S);
}

// ---------------------------------------------------------------------
//  Integer matrix power.  e == 0 -> identity, e < 0 -> inverse^|e|.
// ---------------------------------------------------------------------
static mat mpow(const mat& A, int e) {
    int n = A.n_rows;
    if (e == 0) return eye<mat>(n, n);
    if (e < 0) {
        mat Ai = spd_inv(A), R = Ai;
        for (int i = 1; i < -e; ++i) R = R * Ai;
        return R;
    }
    mat R = A;
    for (int i = 1; i < e; ++i) R = R * A;
    return R;
}

// ---------------------------------------------------------------------
//  "part" matrix and "coeff" scalar for the directional-derivative /
//  equivalence-theorem evaluations (unchanged: operates only on the
//  assembled k x k matrix M and the v x k selection matrix wb).
//
//  pp = 0 :  part = T'(T wb')^{-1} T          coeff = 1
//  pp = 1 :  part = T' T                       coeff = 1/v
//  pp >=2 :  part = T'(T wb')^{p-1} T          coeff = (1/v)^{1/p}
//                                                    * tr(S^p)^{1/p-1}
//  with T = wb * M^{-1}.
// ---------------------------------------------------------------------
static void part_coeff(int pp, const mat& M, const mat& wb,
                       mat& part, double& coeff) {
    mat    temp_inv = wb * spd_inv(M);       // v x k
    double v = (double) wb.n_rows;
    if (pp == 0) {
        part  = temp_inv.t() * spd_inv(temp_inv * wb.t()) * temp_inv;
        coeff = 1.0;
    } else if (pp == 1) {
        part  = temp_inv.t() * temp_inv;
        coeff = 1.0 / v;
    } else {
        mat var_cov = temp_inv * wb.t();    // v x v
        part  = temp_inv.t() * mpow(var_cov, pp - 1) * temp_inv;
        coeff = std::pow(1.0 / v, 1.0 / pp) *
                std::pow(trace(mpow(var_cov, pp)), 1.0 / pp - 1.0);
    }
}

// ---------------------------------------------------------------------
//  Per-point information matrix of candidate column `col` (k x k).
//  Mode 0: outer product of the length-k vector.  Mode 1: reshape the
//  length-k^2 column and symmetrise.
// ---------------------------------------------------------------------
static inline mat point_info(int info_mode, const mat& data, int k, int col) {
    if (info_mode == 0) {
        vec f = data.col(col);
        return f * f.t();
    }
    mat M = reshape(data.col(col), k, k);
    return 0.5 * (M + M.t());
}

// ---------------------------------------------------------------------
//  infor_ind  --  information matrix of a weighted design (candidate
//                 contribution only; the b-scaling is already baked into
//                 `data`).  idx is 1-based.
// ---------------------------------------------------------------------
static mat infor_ind(const std::vector<int>& idx, const vec& w,
                     int info_mode, const mat& data, int k) {
    mat infor(k, k, fill::zeros);
    for (size_t i = 0; i < idx.size(); ++i)
        infor += w(i) * point_info(info_mode, data, k, idx[i] - 1);
    return infor;
}

// ---------------------------------------------------------------------
//  iboss  --  IBOSS subset selection on the covariate-like rows of a
//             k x N information-vector matrix (VECTOR MODE ONLY).
//             X is p x N (the non-intercept rows).  Returns 1-based
//             design-point indices.
// ---------------------------------------------------------------------
static std::vector<int> iboss(const mat& X, int subsize) {
    int p = X.n_rows, N = X.n_cols;
    int cut = (int) std::floor((double) subsize / 2.0 / (double) p);
    std::vector<int> pool(N);
    for (int i = 0; i < N; ++i) pool[i] = i;          // 0-based
    std::vector<int> out;

    for (int row = 0; row < p; ++row) {
        int m = (int) pool.size();
        int c = std::min(cut, m / 2);
        if (c <= 0) break;
        std::vector<int> ord(m);
        for (int i = 0; i < m; ++i) ord[i] = i;
        std::vector<double> vals(m);
        for (int i = 0; i < m; ++i) vals[i] = X(row, pool[i]);
        std::sort(ord.begin(), ord.end(),
                  [&](int a, int b) { return vals[a] < vals[b]; });

        std::vector<int> chosen;
        for (int i = 0; i < c; ++i) chosen.push_back(ord[i]);
        for (int i = 0; i < c; ++i) chosen.push_back(ord[m - c + i]);
        for (int pos : chosen) out.push_back(pool[pos]);

        std::sort(chosen.begin(), chosen.end());
        std::vector<int> next;
        int ci = 0;
        for (int i = 0; i < m; ++i) {
            if (ci < (int) chosen.size() && chosen[ci] == i) { ++ci; continue; }
            next.push_back(pool[i]);
        }
        pool.swap(next);
    }
    int rem = subsize - (int) out.size();
    for (int i = 0; i < rem && i < (int) pool.size(); ++i)
        out.push_back(pool[i]);
    for (int& v : out) v += 1;                        // 1-based
    return out;
}

// ---------------------------------------------------------------------
//  new_weight1_core  --  Newton optimisation of the design weights for a
//                        FIXED support (weight.jl :: new_weight1).
//
//  pinfo holds the per-point k x k information matrices (b-scaling
//  already applied); pinfo.back() is the absorbing point.  infor0 is the
//  existing-design baseline a*I_xi0 (zeros for a single-stage design).
//  This body is representation-agnostic: the only difference between the
//  vector and matrix paths is how pinfo was built.
// ---------------------------------------------------------------------
static void new_weight1_core(int pp, const std::vector<mat>& pinfo,
                             const mat& infor0, const vec& w0, const mat& wb,
                             double r_weight, vec& weight_out, int& stop_out) {
    int n = (int) pinfo.size();
    int v = wb.n_rows;

    const mat& last_infor = pinfo[n - 1];
    std::vector<mat> dinfor(n - 1), infor_blk(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        infor_blk[i] = pinfo[i];
        dinfor[i]    = pinfo[i] - last_infor;
    }

    vec    weight = w0.subvec(0, n - 2);
    double w_diff = 1.0, repli = 1.0, indic = 1.0, delta = 1.0;
    vec    d1w(n - 1, fill::ones);
    mat    d2w(n - 1, n - 1, fill::ones);
    vec    new_weight = weight;

    while (w_diff > 1e-6 && indic > 0.5 && repli < 40.0) {
        mat infor = (r_weight - accu(weight)) * last_infor;
        for (int i = 0; i < n - 1; ++i) infor += weight(i) * infor_blk[i];

        mat temp_inv = spd_inv(infor + infor0);        // k x k
        mat var_cov  = wb * temp_inv * wb.t();          // v x v
        mat inv_var  = mpow(var_cov, pp - 1);           // v x v

        std::vector<mat> part1(n - 1), part2(n - 1);
        std::vector<mat> T1(n - 1), IV2(n - 1);
        for (int i = 0; i < n - 1; ++i) {
            part1[i] =  wb * temp_inv * dinfor[i];                      // v x k
            part2[i] = -wb * temp_inv * dinfor[i] * temp_inv * wb.t();  // v x v
            T1[i]    =  inv_var * part1[i] * temp_inv;                  // v x k
            IV2[i]   =  inv_var * part2[i];                             // v x v
        }
        std::vector<mat> stk;
        std::vector<double> si;
        if (pp > 2) {
            stk.resize(pp - 1);
            stk[0] = eye<mat>(v, v);
            for (int i = 1; i <= pp - 2; ++i) stk[i] = stk[i - 1] * var_cov;
            si.assign(n - 1, 0.0);
            int pm2 = pp - 2, half = pm2 / 2;
            for (int i = 0; i < n - 1; ++i) {
                double s = 0.0;
                if (pm2 == 2 * half) {
                    s += trace(stk[half] * part2[i] * stk[half]);
                    for (int c = 0; c < half; ++c)
                        s += 2.0 * trace(stk[c] * part2[i] * stk[pm2 - c]);
                } else {
                    for (int c = 0; c <= half; ++c)
                        s += 2.0 * trace(stk[c] * part2[i] * stk[pm2 - c]);
                }
                si[i] = s;
            }
        }

        for (int i = 0; i < n - 1; ++i) {
            d1w(i) = trace(IV2[i]);
            for (int j = 0; j <= i; ++j) {
                double val;
                if (pp == 0) {
                    val = 2.0 * accu(T1[i] % part1[j])
                          - accu(IV2[i] % IV2[j].t());
                } else if (pp == 1) {
                    val = 2.0 * accu(T1[i] % part1[j]);
                } else {
                    val = 4.0 * accu(T1[i] % part1[j]) + si[i];
                }
                d2w(i, j) = val;
                d2w(j, i) = val;
            }
        }

        new_weight = weight - delta * (pinv(d2w) * d1w);

        if (new_weight.min() < 0.0 || accu(new_weight) > r_weight) {
            if (delta > 1e-10) delta /= 2.0; else indic = 0.0;
        } else {
            w_diff  = norm(d1w);
            repli  += 1.0;
            weight  = new_weight;
        }
    }

    stop_out = (norm(d1w) > 1e-5) ? 0 : 1;
    weight_out.set_size(n);
    weight_out.subvec(0, n - 2) = new_weight;
    weight_out(n - 1)           = r_weight - accu(new_weight);
}

// Build the per-point information matrices for a 1-based index set.
static std::vector<mat> build_pinfo(const std::vector<int>& ind,
                                    int info_mode, const mat& data, int k) {
    std::vector<mat> pinfo(ind.size());
    for (size_t i = 0; i < ind.size(); ++i)
        pinfo[i] = point_info(info_mode, data, k, ind[i] - 1);
    return pinfo;
}

// ---------------------------------------------------------------------
//  new_weight2  --  optimal weights with zero-weight point removal
//                   (weight.jl :: new_weight2).  Keeps at least min_support
//                   points so the COMBINED information matrix (infor0 + new)
//                   can stay non-singular.  For a single-stage design
//                   min_support = k; for a multistage design the existing
//                   design already supplies rank, so min_support is smaller.
// ---------------------------------------------------------------------
static void new_weight2(int pp, std::vector<int> ind, const mat& infor0,
                        vec w0, const mat& wb, int info_mode, const mat& data,
                        int k, int min_support, std::vector<int>& out_idx, vec& out_w) {
    vec weight;
    if ((int) ind.size() > 1) {
        std::vector<mat> pinfo = build_pinfo(ind, info_mode, data, k);
        int stop;
        new_weight1_core(pp, pinfo, infor0, w0, wb, 1.0, weight, stop);

        while ((int) ind.size() > 1 && weight.min() < 1e-6) {
            if ((int) ind.size() <= min_support) break;
            uword imin = weight.index_min();
            ind.erase(ind.begin() + imin);
            weight.shed_row(imin);
            if ((int) ind.size() > 1) {
                pinfo = build_pinfo(ind, info_mode, data, k);
                int st;
                new_weight1_core(pp, pinfo, infor0, weight, wb, 1.0, weight, st);
            }
        }
    }
    out_idx = ind;
    if ((int) ind.size() == 1) { out_w.set_size(1); out_w(0) = 1.0; }
    else {
        weight  = clamp(weight, 1e-12, datum::inf);
        weight /= accu(weight);
        out_w   = weight;
    }
}

// ---------------------------------------------------------------------
//  directional_deriv  --  the FULL directional-derivative (sensitivity)
//  vector over every candidate column.  opt_infor is the candidate-design
//  contribution (b * I_xi); infor0 is a * I_xi0.  The "part" comes from the
//  COMBINED matrix; the offset trace uses only the candidate contribution
//  (matching the matrix version's offset = tr(K_dyn * I_xi_mat)).  The whole
//  scan is one batched product.  Entry n is d(x_n): at the optimum d(x) <= 0
//  for all candidates, with equality on the support.
// ---------------------------------------------------------------------
static vec directional_deriv(int pp, const mat& opt_infor, const mat& infor0,
                             const mat& wb, int info_mode, const mat& data, int k) {
    mat part; double coeff;
    part_coeff(pp, opt_infor + infor0, wb, part, coeff);
    double diff1 = trace(opt_infor * part);

    rowvec q;
    if (info_mode == 0) {
        q = sum((part * data) % data, 0);             // f_n' part f_n
    } else {
        vec pv = vectorise(part);                     // k^2
        q = (pv.t() * data);                          // <vec(part), vec(I_n)>
    }
    return (q.t() - diff1) * coeff;                   // length-N
}

// ---------------------------------------------------------------------
//  verify_equiv  --  equivalence-theorem check.  Takes the argmax of the
//  directional-derivative vector; max_d <= tol confirms the equivalence
//  theorem.
// ---------------------------------------------------------------------
static void verify_equiv(int pp, const mat& opt_infor, const mat& infor0,
                         const mat& wb, int info_mode, const mat& data, int k,
                         int& out_idx, double& out_sen) {
    vec   dir = directional_deriv(pp, opt_infor, infor0, wb, info_mode, data, k);
    uword idx = dir.index_max();
    out_idx   = (int) idx + 1;
    out_sen   = dir(idx);
}

// ---------------------------------------------------------------------
//  combine  --  merge support points sharing the same index.
// ---------------------------------------------------------------------
static void combine(const std::vector<int>& idx, const vec& w,
                    std::vector<int>& oidx, vec& ow) {
    std::vector<int>    u;
    std::vector<double> uw;
    for (size_t i = 0; i < idx.size(); ++i) {
        auto it = std::find(u.begin(), u.end(), idx[i]);
        if (it == u.end()) { u.push_back(idx[i]); uw.push_back(w(i)); }
        else                 uw[it - u.begin()] += w(i);
    }
    oidx = u;
    ow.set_size(uw.size());
    for (size_t i = 0; i < uw.size(); ++i) ow(i) = uw[i];
}

// ---------------------------------------------------------------------
//  criterion  --  NORMALISED Phi_p criterion value of an approximate
//                 (possibly two-stage) design.
//
//  pp = 0 :  log| Sigma | / v
//  pp > 0 : ( tr( Sigma^p ) / v )^{1/p} ,   Sigma = wb M^{-1} wb',
//  M = infor0 + (b-scaled candidate information).
// ---------------------------------------------------------------------
static double criterion(int pp, const std::vector<int>& idx, const vec& w,
                        int info_mode, const mat& data, int k,
                        const mat& wb, const mat& infor0) {
    mat infor = infor0 + infor_ind(idx, w, info_mode, data, k);
    mat S     = wb * spd_inv(infor) * wb.t();
    double v  = (double) wb.n_rows;
    if (pp == 0) return std::log(det(S)) / v;
    return std::pow(trace(mpow(S, pp)) / v, 1.0 / pp);
}

// =====================================================================
//  criterion_cpp  --  exported wrapper of criterion() for R.
// =====================================================================
// [[Rcpp::export]]
double criterion_cpp(int pp, IntegerVector idx, NumericVector w,
                     int info_mode, const arma::mat& info_data,
                     const arma::mat& wb, const arma::mat& infor0) {
    std::vector<int> id(idx.begin(), idx.end());
    vec ww(w.begin(), w.size());
    int k = (info_mode == 0) ? (int) info_data.n_rows
                             : (int) std::lround(std::sqrt((double) info_data.n_rows));
    return criterion(pp, id, ww, info_mode, info_data, k, wb, infor0);
}

// =====================================================================
//  verify_equiv_cpp  --  directional-derivative scan for a GIVEN design.
//
//  opt_infor is the candidate-design contribution b * I_xi (computed in R
//  from the design's support, which may be off-grid); infor0 = a * I_xi0.
//  Returns the 1-based grid index of the maximiser and the maximum
//  directional derivative (<= tol confirms the equivalence theorem).
// =====================================================================
// [[Rcpp::export]]
List verify_equiv_cpp(int pp, const arma::mat& wb, int info_mode,
                      const arma::mat& info_data, const arma::mat& opt_infor,
                      const arma::mat& infor0) {
    int k = (info_mode == 0) ? (int) info_data.n_rows
                             : (int) std::lround(std::sqrt((double) info_data.n_rows));
    int idx; double sen;
    verify_equiv(pp, opt_infor, infor0, wb, info_mode, info_data, k, idx, sen);
    return List::create(_["index"] = idx, _["max_d"] = sen);
}

// =====================================================================
//  directional_deriv_cpp  --  the FULL directional-derivative (sensitivity)
//  vector over every candidate column, for a GIVEN design.
//
//  Same inputs as verify_equiv_cpp (opt_infor is the b-scaled candidate
//  contribution; infor0 = a * I_xi0), but returns the whole length-N vector
//  instead of just its argmax.  Used by the exact-design construction to pick
//  the most inefficient support point (argmin over support) to remove and the
//  most informative candidate (argmax over all candidates) to add.
// =====================================================================
// [[Rcpp::export]]
arma::vec directional_deriv_cpp(int pp, const arma::mat& wb, int info_mode,
                                const arma::mat& info_data,
                                const arma::mat& opt_infor,
                                const arma::mat& infor0) {
    int k = (info_mode == 0) ? (int) info_data.n_rows
                             : (int) std::lround(std::sqrt((double) info_data.n_rows));
    return directional_deriv(pp, opt_infor, infor0, wb, info_mode, info_data, k);
}

// =====================================================================
//  optimize_weights_cpp  --  re-optimise the weights on a FIXED support
//  (used by the optional neighbourhood-merge step, where the merged
//  support points are off-grid and supplied as their own info_data).
//  Returns the surviving 1-based indices, their weights, and the
//  normalised criterion value.
// =====================================================================
// [[Rcpp::export]]
List optimize_weights_cpp(int pp, const arma::mat& wb, int info_mode,
                          const arma::mat& info_data, const arma::mat& infor0,
                          int min_support) {
    int k = (info_mode == 0) ? (int) info_data.n_rows
                             : (int) std::lround(std::sqrt((double) info_data.n_rows));
    int m = (int) info_data.n_cols;
    std::vector<int> ind(m);
    for (int i = 0; i < m; ++i) ind[i] = i + 1;          // 1-based, all points
    vec w0(m); w0.fill(1.0 / (double) m);

    std::vector<int> sup; vec sw;
    new_weight2(pp, ind, infor0, w0, wb, info_mode, info_data, k, min_support, sup, sw);
    std::vector<int> cidx; vec cw;
    combine(sup, sw, cidx, cw);

    IntegerVector r_idx(cidx.begin(), cidx.end());
    NumericVector r_w(cw.begin(), cw.end());
    return List::create(
        _["index"]  = r_idx,
        _["weight"] = r_w,
        _["value"]  = criterion(pp, cidx, cw, info_mode, info_data, k, wb, infor0));
}

// =====================================================================
//  appro_opt_cpp  --  the unified full solver.
//
//  pp        criterion (0 = D, 1 = A, >=2 = Phi_p).
//  wb        v x k selection matrix dg/dtheta.
//  info_mode 0 = vector (info_data k x N), 1 = matrix (info_data k^2 x N).
//  info_data candidate information, ALREADY b-scaled by the R layer.
//  infor0    a * I_xi0 (k x k); zeros for a single-stage design.
//  init_idx  1-based initial support; if empty, IBOSS is run (vector mode
//            only -- the R layer always supplies init_idx for matrix mode).
//
//  Returns: index (1-based support), weight, sensitivity (max directional
//  derivative), iter (outer iterations), value (normalised Phi_p).
// =====================================================================
// [[Rcpp::export]]
List appro_opt_cpp(int pp, const arma::mat& wb, int info_mode,
                   const arma::mat& info_data, const arma::mat& infor0,
                   IntegerVector init_idx, int min_support,
                   int max_iter, double tol, bool verbose) {
    int k = (info_mode == 0) ? (int) info_data.n_rows
                             : (int) std::lround(std::sqrt((double) info_data.n_rows));

    // ---- 1. initial support ----------------------------------------
    std::vector<int> idx0;
    if (init_idx.size() > 0) {
        idx0.assign(init_idx.begin(), init_idx.end());        // 1-based
    } else {
        if (info_mode != 0)
            stop("matrix mode requires an explicit init_idx");
        idx0 = iboss(info_data.rows(1, k - 1), 2 * (k - 1));
    }

    // ---- 2. optimal weights on the initial support -----------------
    vec w0(idx0.size());
    w0.fill(1.0 / idx0.size());
    std::vector<int> sup; vec sw;
    new_weight2(pp, idx0, infor0, w0, wb, info_mode, info_data, k, min_support, sup, sw);

    std::vector<int> cidx; vec cw;
    combine(sup, sw, cidx, cw);

    mat opt_infor = infor_ind(cidx, cw, info_mode, info_data, k);
    int    vidx; double sen;
    verify_equiv(pp, opt_infor, infor0, wb, info_mode, info_data, k, vidx, sen);

    int iter = 1;
    // ---- 3-4. add the best point, re-optimise, until d_p <= tol -----
    while (sen > tol && iter < max_iter) {
        std::vector<int> newx = cidx;
        newx.push_back(vidx);
        vec new_w0(cw.n_elem + 1);
        new_w0.subvec(0, cw.n_elem - 1) = cw;
        new_w0(cw.n_elem) = 0.0;

        new_weight2(pp, newx, infor0, new_w0, wb, info_mode, info_data, k, min_support, sup, sw);
        combine(sup, sw, cidx, cw);
        opt_infor = infor_ind(cidx, cw, info_mode, info_data, k);
        verify_equiv(pp, opt_infor, infor0, wb, info_mode, info_data, k, vidx, sen);
        ++iter;
        if (verbose)
            Rcout << "  iter " << iter << "   support " << cidx.size()
                  << "   max d = " << sen << "\n";
    }

    // ---- sort support by index for a tidy result -------------------
    std::vector<int> ord(cidx.size());
    for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](int a, int b) { return cidx[a] < cidx[b]; });
    IntegerVector r_idx(cidx.size());
    NumericVector r_w(cidx.size());
    std::vector<int> s_idx(cidx.size());
    vec              s_w(cidx.size());
    for (size_t i = 0; i < ord.size(); ++i) {
        r_idx[i] = cidx[ord[i]];  s_idx[i] = cidx[ord[i]];
        r_w[i]   = cw(ord[i]);    s_w(i)   = cw(ord[i]);
    }

    return List::create(
        _["index"]       = r_idx,
        _["weight"]      = r_w,
        _["sensitivity"] = sen,
        _["iter"]        = iter,
        _["value"]       = criterion(pp, s_idx, s_w, info_mode, info_data, k,
                                     wb, infor0));
}
