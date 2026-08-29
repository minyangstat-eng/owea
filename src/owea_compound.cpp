// =====================================================================
//  owea_compound.cpp -- C++ core for the COMPOUND (multi-model)
//                       criterion.
//
//  This file is ADDITIVE: it defines its own helpers (all `static`, all
//  prefixed `cmp_`) and its own exported entry points (`compound_*`).
//  It does not touch owea_engine.cpp and does not change the behaviour
//  of any existing function.
//
//  THE CRITERION
//  -------------
//  J components share ONE design xi.  Component j has its own model, so
//  its own per-point information I_j(x), its own parameter dimension
//  k_j, its own existing-design information M_{0,j}, its own quantity of
//  interest G_j (v_j x k_j) and its own type (D or A):
//
//      M_j(xi) = a M_{0,j} + b M_{xi,j},   S_j = G_j M_j^{-1} G_j'
//
//      Psi_j = det(S_j^{-1})^{1/v_j}   (type D, pp = 0)
//      Psi_j = v_j / tr(S_j)           (type A, pp = 1)
//
//      Psi_alpha = sum_j at_j Psi_j ,  at_j = alpha_j / Psi_j^*  (or
//                                      alpha_j for the raw variant).
//
//  THE ONE OBJECT EVERYTHING RUNS THROUGH
//  --------------------------------------
//  Because the components live in different spaces the single "part"
//  matrix of the one-model engine does not exist.  It is replaced by a
//  LIST of k_j x k_j matrices
//
//      K_j = at_j * c_j * P_j ,
//      P_j = T_j' S_j^{-1} T_j (D)   or   T_j' T_j (A),   T_j = G_j M_j^{-1}
//      c_j = Psi_j / v_j       (D)   or   Psi_j^2 / v_j   (A)
//
//  with the Euler identity  sum_j tr(K_j M_j) = Psi_alpha  (each Psi_j is
//  positively homogeneous of degree one in M_j).  Then
//
//      sensitivity   psi(x)   = sum_j [ tr(K_j I_j(x)) - tr(K_j M_{xi,j}) ]
//      gradient      d/du_i   = sum_j tr(K_j D_{j,i}) ,  D_{j,i}=I_j(x_i)-I_j(x_m)
//      Hessian       d2/du_i du_l = sum_j at_j H_{j,il}  (see cmp_grad_hess)
//
//  SCALING CONVENTION (identical to owea_engine.cpp)
//  ------------------------------------------------
//  The R layer PRE-SCALES each component's candidate information by b
//  (vectors by sqrt(b), matrices by b) and passes infor0_j = a M_{0,j}.
//  So every per-point matrix built here already carries b, the assembled
//  `infor0_j + infor_ind_j` equals M_j exactly, and the b factors in the
//  sensitivity/gradient/Hessian above are absorbed automatically.
//
//  Conventions
//   * design-point indices are 1-based across the R boundary;
//   * info_mode 0 = information VECTOR (k x N, I = f f'),
//     info_mode 1 = information MATRIX (k*k x N, column = vec(I));
//   * pp = 0 -> D, pp = 1 -> A (no other value is accepted).
// =====================================================================

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------
//  cmp_spd_inv -- robust inverse of a symmetric positive-(semi)definite
//  matrix: Cholesky, then an increasing ridge, then pinv.  Mirrors the
//  behaviour of spd_inv() in owea_engine.cpp but is a separate symbol.
// ---------------------------------------------------------------------
static mat cmp_spd_inv(const mat& A) {
    mat B = 0.5 * (A + A.t());
    mat R, out;
    if (chol(R, B)) {
        mat Ri = inv(trimatu(R));
        return Ri * Ri.t();
    }
    double ridge = 1e-12 * (std::abs(trace(B)) / std::max(1u, (unsigned) B.n_rows) + 1.0);
    for (int t = 0; t < 12; ++t) {
        mat Bt = B + ridge * eye<mat>(B.n_rows, B.n_cols);
        if (chol(R, Bt)) {
            mat Ri = inv(trimatu(R));
            return Ri * Ri.t();
        }
        ridge *= 10.0;
    }
    if (pinv(out, B)) return out;
    return eye<mat>(B.n_rows, B.n_cols);
}

// log|S| via Cholesky; ok = false when S is not positive definite.
static double cmp_logdet(const mat& S, bool& ok) {
    mat R;
    ok = chol(R, 0.5 * (S + S.t()));
    if (!ok) return 0.0;
    return 2.0 * accu(log(R.diag()));
}

// ---------------------------------------------------------------------
//  One component: its model data and its criterion.
// ---------------------------------------------------------------------
struct CmpComp {
    int    pp;         // 0 = D, 1 = A
    int    info_mode;  // 0 = vector, 1 = matrix
    int    k;          // parameter dimension of THIS model
    int    v;          // rows of wb
    double at;         // tilde alpha_j
    mat    data;       // k x N  or  k*k x N, already b-scaled
    mat    wb;         // v x k
    mat    infor0;     // a * M_{0,j}
};

// Per-component quantities at the current design.
struct CmpState {
    mat    M, Mi, T, S, Si, P, K, opt_infor;
    double psi;   // Psi_j
    double c;     // c_j
    bool   ok;    // M_j and S_j usable
};

// ---------------------------------------------------------------------
//  Per-point information matrix of candidate column `col`, model j.
// ---------------------------------------------------------------------
static inline mat cmp_point_info(const CmpComp& C, int col) {
    if (C.info_mode == 0) {
        vec f = C.data.col(col);
        return f * f.t();
    }
    mat M = reshape(C.data.col(col), C.k, C.k);
    return 0.5 * (M + M.t());
}

// Candidate contribution sum_i w_i I_j(x_i) (b already inside `data`).
// idx is 1-based.
static mat cmp_infor_ind(const CmpComp& C, const std::vector<int>& idx,
                         const vec& w) {
    mat out(C.k, C.k, fill::zeros);
    for (size_t i = 0; i < idx.size(); ++i)
        out += w(i) * cmp_point_info(C, idx[i] - 1);
    return out;
}

// ---------------------------------------------------------------------
//  cmp_state -- assemble M_j = infor0_j + opt_infor_j and everything
//  derived from it, including K_j = at_j c_j P_j.
// ---------------------------------------------------------------------
static CmpState cmp_state(const CmpComp& C, const mat& opt_infor) {
    CmpState S;
    S.opt_infor = opt_infor;
    S.M  = C.infor0 + opt_infor;
    S.Mi = cmp_spd_inv(S.M);
    S.T  = C.wb * S.Mi;                 // v x k
    S.S  = S.T * C.wb.t();              // v x v
    S.ok = true;

    if (C.pp == 0) {                    // ---- D ----
        bool ok = false;
        double ld = cmp_logdet(S.S, ok);
        if (!ok) { S.ok = false; ld = 0.0; }
        S.psi = std::exp(-ld / (double) C.v);
        S.Si  = cmp_spd_inv(S.S);
        S.P   = S.T.t() * S.Si * S.T;
        S.c   = S.psi / (double) C.v;
    } else {                            // ---- A ----
        double trS = trace(S.S);
        if (!(trS > 0.0) || !std::isfinite(trS)) { S.ok = false; trS = 1.0; }
        S.psi = (double) C.v / trS;
        S.Si  = eye<mat>(C.v, C.v);     // unused for A
        S.P   = S.T.t() * S.T;
        S.c   = S.psi * S.psi / (double) C.v;
    }
    if (!std::isfinite(S.psi)) { S.ok = false; S.psi = 0.0; }
    S.K = C.at * S.c * S.P;
    return S;
}

static std::vector<CmpState> cmp_states(const std::vector<CmpComp>& C,
                                        const std::vector<int>& idx,
                                        const vec& w) {
    std::vector<CmpState> St(C.size());
    for (size_t j = 0; j < C.size(); ++j)
        St[j] = cmp_state(C[j], cmp_infor_ind(C[j], idx, w));
    return St;
}

// Psi_alpha = sum_j at_j Psi_j.
static double cmp_criterion(const std::vector<CmpComp>& C,
                            const std::vector<CmpState>& St) {
    double out = 0.0;
    for (size_t j = 0; j < C.size(); ++j) out += C[j].at * St[j].psi;
    return out;
}

// ---------------------------------------------------------------------
//  cmp_dirderiv -- sensitivity over EVERY candidate:
//      psi(x_n) = sum_j [ tr(K_j I_j(x_n)) - tr(K_j M_{xi,j}) ]
//  One batched product per component; the offset is subtracted once.
//  At the optimum psi(x) <= 0 for all n, with equality on the support.
// ---------------------------------------------------------------------
static vec cmp_dirderiv(const std::vector<CmpComp>& C,
                        const std::vector<CmpState>& St) {
    const uword N = C[0].data.n_cols;
    vec out(N, fill::zeros);
    double offset = 0.0;
    for (size_t j = 0; j < C.size(); ++j) {
        rowvec q;
        if (C[j].info_mode == 0) {
            q = sum((St[j].K * C[j].data) % C[j].data, 0);   // f' K f
        } else {
            vec pv = vectorise(St[j].K);
            q = pv.t() * C[j].data;                          // <vec(K), vec(I)>
        }
        out    += q.t();
        offset += trace(St[j].K * St[j].opt_infor);
    }
    return out - offset;
}

static void cmp_verify(const std::vector<CmpComp>& C,
                       const std::vector<CmpState>& St,
                       int& out_idx, double& out_sen) {
    vec  d   = cmp_dirderiv(C, St);
    uword im = d.index_max();
    out_idx  = (int) im + 1;
    out_sen  = d(im);
}

// ---------------------------------------------------------------------
//  cmp_grad_hess -- gradient and Hessian of Psi_alpha with respect to the
//  free weights u = (w_1, ..., w_{n-1}),  w_n = r_weight - sum_i u_i.
//
//      dPsi/du_i        = sum_j tr(K_j D_{j,i})
//      d2Psi/du_i du_l  = sum_j at_j H_{j,il} ,
//
//      D-type: H = (Psi/v) [ tr(P D_i) tr(P D_l)/v
//                            - 2 tr(P D_i M^{-1} D_l)
//                            + tr(S^{-1} A_i S^{-1} A_l) ]
//      A-type: H = (2 Psi^2/v) [ (Psi/v) tr(P D_i) tr(P D_l)
//                                - tr(P D_i M^{-1} D_l) ]
//
//  with A_i = T D_i T'.  All D_{j,i} carry b, so no b appears here.
// ---------------------------------------------------------------------
static void cmp_grad_hess(const std::vector<CmpComp>& C,
                          const std::vector<CmpState>& St,
                          const std::vector< std::vector<mat> >& D,  // D[j][i]
                          vec& g, mat& H) {
    const int J  = (int) C.size();
    const int mm = (int) D[0].size();          // n - 1
    g.zeros(mm);
    H.zeros(mm, mm);

    for (int j = 0; j < J; ++j) {
        const CmpState& s = St[j];
        const double v    = (double) C[j].v;

        std::vector<mat>    E(mm);             // E_i = P D_i M^{-1}
        std::vector<mat>    SA(mm);            // SA_i = S^{-1} A_i   (D only)
        std::vector<double> tr_i(mm);          // tr(P D_i)

        for (int i = 0; i < mm; ++i) {
            mat PD  = s.P * D[j][i];
            tr_i[i] = trace(PD);
            E[i]    = PD * s.Mi;
            if (C[j].pp == 0) {
                mat A = s.T * D[j][i] * s.T.t();
                SA[i] = s.Si * A;
            }
            g(i) += C[j].at * s.c * tr_i[i];   // = tr(K_j D_{j,i})
        }

        const double pre = (C[j].pp == 0) ? (s.psi / v)
                                          : (2.0 * s.psi * s.psi / v);
        for (int i = 0; i < mm; ++i) {
            for (int l = 0; l <= i; ++l) {
                double mid = accu(E[i] % D[j][l]);          // tr(P D_i M^-1 D_l)
                double val;
                if (C[j].pp == 0) {
                    double quad = accu(SA[i] % SA[l].t());  // tr(S^-1A_i S^-1A_l)
                    val = pre * (tr_i[i] * tr_i[l] / v - 2.0 * mid + quad);
                } else {
                    val = pre * ((s.psi / v) * tr_i[i] * tr_i[l] - mid);
                }
                val *= C[j].at;
                H(i, l) += val;
                if (l != i) H(l, i) += val;
            }
        }
    }
}

// ---------------------------------------------------------------------
//  cmp_weights_core -- damped Newton on the weights of a FIXED support.
//
//      u <- u - delta * H^{-1} g          (H = Hessian of Psi_alpha)
//
//  H is negative semi-definite (Psi_alpha is concave), so -H^{-1}g is an
//  ascent direction and this is the maximisation step.  delta is halved
//  until the new weights stay in the simplex, every M_j stays usable, and
//  Psi_alpha does not decrease.
// ---------------------------------------------------------------------
static void cmp_weights_core(const std::vector<CmpComp>& C,
                             const std::vector< std::vector<mat> >& pinfo, // [j][i]
                             const vec& w0, double r_weight,
                             vec& weight_out, int& stop_out) {
    const int J = (int) C.size();
    const int n = (int) pinfo[0].size();

    if (n < 2) {
        weight_out.set_size(std::max(n, 1));
        weight_out.fill(r_weight);
        stop_out = 1;
        return;
    }
    const int mm = n - 1;

    // D[j][i] = I_j(x_i) - I_j(x_n)
    std::vector< std::vector<mat> > D(J, std::vector<mat>(mm));
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < mm; ++i)
            D[j][i] = pinfo[j][i] - pinfo[j][n - 1];

    // assemble the states for a given free-weight vector
    auto states_of = [&](const vec& u) {
        std::vector<CmpState> St(J);
        double wlast = r_weight - accu(u);
        for (int j = 0; j < J; ++j) {
            mat oi = wlast * pinfo[j][n - 1];
            for (int i = 0; i < mm; ++i) oi += u(i) * pinfo[j][i];
            St[j] = cmp_state(C[j], oi);
        }
        return St;
    };
    auto usable = [&](const std::vector<CmpState>& St) {
        for (int j = 0; j < J; ++j) if (!St[j].ok) return false;
        return true;
    };

    vec u = w0.subvec(0, mm - 1);
    std::vector<CmpState> St = states_of(u);
    double cur = cmp_criterion(C, St);

    vec g; mat H;
    double gnorm = datum::inf;

    for (int it = 0; it < 200; ++it) {
        cmp_grad_hess(C, St, D, g, H);
        gnorm = norm(g);
        if (gnorm < 1e-12) break;

        // Two candidate ascent directions, tried in order:
        //   (1) the Newton direction -H^{-1} g  (H is negative semi-definite,
        //       so this ascends), and
        //   (2) plain gradient ascent, as a fallback when H is so
        //       ill-conditioned that every damped Newton trial leaves the
        //       simplex.  Without (2) the iteration can stall with a non-zero
        //       gradient and report a spuriously unconverged design.
        std::vector<vec> dirs;
        vec step;
        if (solve(step, H, g, solve_opts::no_approx) && step.is_finite())
            dirs.push_back(-step);                 // u + delta * (-H^{-1} g)
        else {
            vec ps = pinv(H) * g;
            if (ps.is_finite()) dirs.push_back(-ps);
        }
        {   // scale the gradient direction to a sensible first trial length
            double gn = norm(g);
            if (gn > 0.0 && std::isfinite(gn)) dirs.push_back(g / gn);
        }

        bool moved = false;
        for (size_t d = 0; d < dirs.size() && !moved; ++d) {
            double delta = 1.0;
            for (int t = 0; t < 60; ++t) {
                vec un = u + delta * dirs[d];
                if (un.min() >= 0.0 && accu(un) <= r_weight) {
                    std::vector<CmpState> Sn = states_of(un);
                    if (usable(Sn)) {
                        double val = cmp_criterion(C, Sn);
                        if (val > cur) {
                            u = un; St = Sn; cur = val; moved = true;
                            break;
                        }
                    }
                }
                delta *= 0.5;
            }
        }
        if (!moved) break;
    }

    weight_out.set_size(n);
    weight_out.subvec(0, mm - 1) = u;
    weight_out(n - 1) = r_weight - accu(u);
    stop_out = (gnorm > 1e-5) ? 0 : 1;
}

// Build the per-component per-point information for a 1-based index set.
static std::vector< std::vector<mat> >
cmp_build_pinfo(const std::vector<CmpComp>& C, const std::vector<int>& idx) {
    std::vector< std::vector<mat> > pinfo(C.size(), std::vector<mat>(idx.size()));
    for (size_t j = 0; j < C.size(); ++j)
        for (size_t i = 0; i < idx.size(); ++i)
            pinfo[j][i] = cmp_point_info(C[j], idx[i] - 1);
    return pinfo;
}

// ---------------------------------------------------------------------
//  cmp_weights2 -- optimal weights with zero-weight point removal.
//  Dropping points that reach the boundary is essential: with an active
//  constraint the Newton step keeps pointing out of the simplex and the
//  damping drives delta to zero with a non-zero gradient.
//  Keeps at least min_support points.
// ---------------------------------------------------------------------
static void cmp_weights2(const std::vector<CmpComp>& C,
                         std::vector<int> ind, const vec& w0,
                         int min_support,
                         std::vector<int>& out_idx, vec& out_w) {
    vec weight;
    if (ind.size() > 1) {
        std::vector< std::vector<mat> > pinfo = cmp_build_pinfo(C, ind);
        int stop;
        cmp_weights_core(C, pinfo, w0, 1.0, weight, stop);

        while (ind.size() > 1 && weight.min() < 1e-6) {
            if ((int) ind.size() <= min_support) break;
            uword imin = weight.index_min();
            ind.erase(ind.begin() + imin);
            weight.shed_row(imin);
            if (ind.size() > 1) {
                weight /= accu(weight);
                pinfo = cmp_build_pinfo(C, ind);
                int st;
                cmp_weights_core(C, pinfo, weight, 1.0, weight, st);
            }
        }
    }
    out_idx = ind;
    if (ind.size() == 1) { out_w.set_size(1); out_w(0) = 1.0; }
    else {
        weight  = clamp(weight, 1e-12, datum::inf);
        weight /= accu(weight);
        out_w   = weight;
    }
}

// Merge duplicated support indices (same helper role as combine()).
static void cmp_combine(const std::vector<int>& idx, const vec& w,
                        std::vector<int>& oidx, vec& ow) {
    std::vector<int>    u;
    std::vector<double> uw;
    for (size_t i = 0; i < idx.size(); ++i) {
        std::vector<int>::iterator it = std::find(u.begin(), u.end(), idx[i]);
        if (it == u.end()) { u.push_back(idx[i]); uw.push_back(w(i)); }
        else                 uw[it - u.begin()] += w(i);
    }
    oidx = u;
    ow.set_size(uw.size());
    for (size_t i = 0; i < uw.size(); ++i) ow(i) = uw[i];
}

// ---------------------------------------------------------------------
//  Unpack the R-side lists into components.
// ---------------------------------------------------------------------
static std::vector<CmpComp> cmp_unpack(List info_data, IntegerVector info_mode,
                                       List wb, List infor0, IntegerVector pp,
                                       NumericVector at) {
    const int J = info_data.size();
    if (wb.size() != J || infor0.size() != J || pp.size() != J ||
        at.size() != J || info_mode.size() != J)
        stop("compound: component lists must all have the same length.");

    std::vector<CmpComp> C(J);
    uword N = 0;
    for (int j = 0; j < J; ++j) {
        CmpComp& c = C[j];
        c.pp        = pp[j];
        c.info_mode = info_mode[j];
        c.at        = at[j];
        c.data      = as<mat>(info_data[j]);
        c.wb        = as<mat>(wb[j]);
        c.infor0    = as<mat>(infor0[j]);
        if (c.pp != 0 && c.pp != 1)
            stop("compound: each component's p must be 0 (D) or 1 (A).");
        c.k = (c.info_mode == 0)
                ? (int) c.data.n_rows
                : (int) std::lround(std::sqrt((double) c.data.n_rows));
        c.v = (int) c.wb.n_rows;
        if ((int) c.wb.n_cols != c.k)
            stop("compound: wb for a component has the wrong number of columns.");
        if ((int) c.infor0.n_rows != c.k || (int) c.infor0.n_cols != c.k)
            stop("compound: infor0 for a component has the wrong dimension.");
        if (j == 0) N = c.data.n_cols;
        else if (c.data.n_cols != N)
            stop("compound: all components must use the same candidate set.");
    }
    return C;
}

// =====================================================================
//  EXPORTED ENTRY POINTS
// =====================================================================

//' Compound criterion value (internal).
//' @keywords internal
// [[Rcpp::export]]
double compound_criterion_cpp(List info_data, IntegerVector info_mode, List wb,
                              List infor0, IntegerVector pp, NumericVector at,
                              IntegerVector idx, NumericVector w) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    std::vector<CmpState> St = cmp_states(C, id, ww);
    return cmp_criterion(C, St);
}

//' Compound per-component criterion values (internal).
//' @keywords internal
// [[Rcpp::export]]
NumericVector compound_psi_cpp(List info_data, IntegerVector info_mode, List wb,
                               List infor0, IntegerVector pp, NumericVector at,
                               IntegerVector idx, NumericVector w) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    std::vector<CmpState> St = cmp_states(C, id, ww);
    NumericVector out(C.size());
    for (size_t j = 0; j < C.size(); ++j) out[j] = St[j].psi;
    return out;
}

//' Compound sensitivity over every candidate (internal).
//' @keywords internal
// [[Rcpp::export]]
NumericVector compound_dirderiv_cpp(List info_data, IntegerVector info_mode,
                                    List wb, List infor0, IntegerVector pp,
                                    NumericVector at, IntegerVector idx,
                                    NumericVector w) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    std::vector<CmpState> St = cmp_states(C, id, ww);
    vec d = cmp_dirderiv(C, St);
    return NumericVector(d.begin(), d.end());
}

//' Compound equivalence-theorem check (internal).
//' @keywords internal
// [[Rcpp::export]]
List compound_verify_cpp(List info_data, IntegerVector info_mode, List wb,
                         List infor0, IntegerVector pp, NumericVector at,
                         IntegerVector idx, NumericVector w) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    std::vector<CmpState> St = cmp_states(C, id, ww);
    int vi; double sen;
    cmp_verify(C, St, vi, sen);

    // the Euler identity sum_j tr(K_j M_j) = Psi_alpha, for diagnostics
    double euler = 0.0;
    for (size_t j = 0; j < C.size(); ++j) euler += trace(St[j].K * St[j].M);

    return List::create(_["index"]  = vi,
                        _["max_d"]  = sen,
                        _["value"]  = cmp_criterion(C, St),
                        _["euler"]  = euler);
}

//' Compound gradient and Hessian in the free weights (internal).
//' @keywords internal
// [[Rcpp::export]]
List compound_grad_hess_cpp(List info_data, IntegerVector info_mode, List wb,
                            List infor0, IntegerVector pp, NumericVector at,
                            IntegerVector idx, NumericVector w) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    const int n = (int) id.size();
    if (n < 2) stop("compound: need at least two support points.");

    std::vector< std::vector<mat> > pinfo = cmp_build_pinfo(C, id);
    const int J = (int) C.size(), mm = n - 1;
    std::vector< std::vector<mat> > D(J, std::vector<mat>(mm));
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < mm; ++i) D[j][i] = pinfo[j][i] - pinfo[j][n - 1];

    std::vector<CmpState> St = cmp_states(C, id, ww);
    vec g; mat H;
    cmp_grad_hess(C, St, D, g, H);
    return List::create(_["gradient"] = NumericVector(g.begin(), g.end()),
                        _["hessian"]  = wrap(H));
}

//' Optimal compound weights on a fixed support (internal).
//' @keywords internal
// [[Rcpp::export]]
List compound_weights_cpp(List info_data, IntegerVector info_mode, List wb,
                          List infor0, IntegerVector pp, NumericVector at,
                          IntegerVector idx, NumericVector w, int min_support) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    std::vector<int> id(idx.begin(), idx.end());
    vec ww = as<vec>(w);
    std::vector<int> oi; vec ow;
    cmp_weights2(C, id, ww, min_support, oi, ow);
    std::vector<CmpState> St = cmp_states(C, oi, ow);
    return List::create(_["index"]  = IntegerVector(oi.begin(), oi.end()),
                        _["weight"] = NumericVector(ow.begin(), ow.end()),
                        _["value"]  = cmp_criterion(C, St));
}

//' Compound optimal design on a fixed candidate set (internal).
//'
//' Exchange loop mirroring \code{appro_opt_cpp}: optimise the weights of the
//' current support, add the candidate of maximum compound sensitivity, repeat
//' until the sensitivity falls to \code{tol}.
//' @keywords internal
// [[Rcpp::export]]
List compound_appro_opt_cpp(List info_data, IntegerVector info_mode, List wb,
                            List infor0, IntegerVector pp, NumericVector at,
                            IntegerVector init_idx, int min_support,
                            int max_iter, double tol, bool verbose) {
    std::vector<CmpComp> C = cmp_unpack(info_data, info_mode, wb, infor0, pp, at);
    if (init_idx.size() == 0)
        stop("compound: an explicit init_idx is required.");

    // ---- 1. initial support -----------------------------------------
    std::vector<int> idx0(init_idx.begin(), init_idx.end());
    vec w0(idx0.size());
    w0.fill(1.0 / idx0.size());

    // ---- 2. optimal weights on it -----------------------------------
    std::vector<int> sup; vec sw;
    cmp_weights2(C, idx0, w0, min_support, sup, sw);

    std::vector<int> cidx; vec cw;
    cmp_combine(sup, sw, cidx, cw);

    std::vector<CmpState> St = cmp_states(C, cidx, cw);
    int vidx; double sen;
    cmp_verify(C, St, vidx, sen);

    int iter = 1;
    // ---- 3-4. add the best point, re-optimise, until psi <= tol ------
    while (sen > tol && iter < max_iter) {
        std::vector<int> newx = cidx;
        vec new_w0;
        bool present = std::find(cidx.begin(), cidx.end(), vidx) != cidx.end();
        if (present) {
            // already in the support: the weights are simply not optimal yet;
            // re-optimising from the current point is the right move.
            new_w0 = cw;
        } else {
            newx.push_back(vidx);
            new_w0.set_size(cw.n_elem + 1);
            new_w0.subvec(0, cw.n_elem - 1) = cw * (1.0 - 1e-3);
            new_w0(cw.n_elem) = 1e-3;
        }

        cmp_weights2(C, newx, new_w0, min_support, sup, sw);
        cmp_combine(sup, sw, cidx, cw);
        St = cmp_states(C, cidx, cw);
        double prev = sen;
        cmp_verify(C, St, vidx, sen);
        ++iter;
        if (verbose)
            Rcout << "  iter " << iter << "   support " << cidx.size()
                  << "   max d = " << sen << "\n";
        if (present && sen >= prev - 1e-15) break;   // no further progress
    }

    // ---- sort the support by index ----------------------------------
    std::vector<int> ord(cidx.size());
    for (size_t i = 0; i < ord.size(); ++i) ord[i] = (int) i;
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

    std::vector<CmpState> Sf = cmp_states(C, s_idx, s_w);
    NumericVector psi_j(C.size());
    for (size_t j = 0; j < C.size(); ++j) psi_j[j] = Sf[j].psi;

    return List::create(_["index"]       = r_idx,
                        _["weight"]      = r_w,
                        _["sensitivity"] = sen,
                        _["iter"]        = iter,
                        _["value"]       = cmp_criterion(C, Sf),
                        _["psi"]         = psi_j);
}
