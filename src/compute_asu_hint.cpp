// src/compute_asu_hint.cpp
// Warm-start hint algorithms for the CP-SAT ASU solver, ported from
// inst/python/asu_cpsat.py (greedy_snake_hint, _articulation_points,
// reverse_prune_hint).  All indices are 0-based throughout, matching
// Python's convention.  The exported functions are designed to be called
// from R (via reticulate or directly) before launching the Python solver so
// that the pre-computed hint can be passed as an initial solution.
//
// Exported R functions:
//   cpp_articulation_points(nb_local, selected)   -> integer vector (0-based cut vertices)
//   cpp_greedy_snake_hint(nb_local, u_g, E_g, P_g, tau, pop_thresh, root_local)
//   cpp_reverse_prune_hint(nb_local, u_g, E_g, P_g, tau, pop_thresh, root_local)
//   cpp_compute_asu_hint(nb_local, u_g, E_g, P_g, tau, pop_thresh, root_local)
//
// All nb_local entries are 0-based integer vectors.

#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// -----------------------------------------------------------------------
// Internal helpers
// -----------------------------------------------------------------------

// Parse the neighbor List into a flat C++ structure.
static void parse_nb(const List& nb_local, std::vector<std::vector<int>>& nb) {
    int N = nb_local.size();
    nb.resize(N);
    for (int i = 0; i < N; ++i) {
        IntegerVector row = nb_local[i];
        nb[i].assign(row.begin(), row.end());
    }
}

static inline double ur_of(double u, double e) {
    return (u + e > 0.0) ? u / (u + e) : 0.0;
}

// -----------------------------------------------------------------------
// Iterative Tarjan articulation-point finder (mirrors _articulation_points)
// Restricted to the induced subgraph on nodes where selected[i] == true.
// Returns 0-based indices of cut vertices.
// -----------------------------------------------------------------------
static std::unordered_set<int> art_pts_impl(
    const std::vector<std::vector<int>>& nb,
    const std::vector<bool>& selected
) {
    int N = static_cast<int>(nb.size());
    std::vector<int>  disc(N, -1), low(N, 0), par(N, -1), pos(N, 0), root_ch(N, 0);
    std::vector<bool> skip_par(N, false);
    std::unordered_set<int> art;
    int timer = 0;

    for (int start = 0; start < N; ++start) {
        if (!selected[start] || disc[start] != -1) continue;
        disc[start] = low[start] = timer++;
        std::vector<int> stk = {start};

        while (!stk.empty()) {
            int u = stk.back();
            bool pushed = false;
            while (pos[u] < static_cast<int>(nb[u].size())) {
                int w = nb[u][pos[u]++];
                if (!selected[w]) continue;
                if (w == par[u] && !skip_par[u]) { skip_par[u] = true; continue; }
                if (disc[w] == -1) {
                    par[w] = u;
                    if (par[u] == -1) ++root_ch[start];  // u is root, w is root's direct child
                    disc[w] = low[w] = timer++;
                    stk.push_back(w);
                    pushed = true;
                    break;
                }
                low[u] = std::min(low[u], disc[w]);
            }
            if (!pushed) {
                stk.pop_back();
                if (!stk.empty()) {
                    int p = stk.back();
                    low[p] = std::min(low[p], low[u]);
                    if (par[p] != -1 && low[u] >= disc[p]) art.insert(p);
                }
            }
        }
        if (root_ch[start] > 1) art.insert(start);
    }
    return art;
}

// -----------------------------------------------------------------------
// Greedy-snake single-group expand from one seed (mirrors _expand closure).
// remaining[] is read-only here; caller marks assigned nodes False after.
// -----------------------------------------------------------------------
static std::vector<int> expand_impl(
    int seed,
    const std::vector<std::vector<int>>& nb,
    const std::vector<double>& u,
    const std::vector<double>& E,
    const std::vector<double>& P,
    double tau,
    int pop_thresh,
    const std::vector<bool>& remaining   // NOT modified inside
) {
    std::unordered_set<int> sel_set = {seed};
    double su = u[seed], se = E[seed], sp = P[seed];

    std::unordered_set<int> frontier;
    for (int w : nb[seed]) {
        if (remaining[w]) frontier.insert(w);
    }

    while (!frontier.empty()) {
        int best = -1;
        double best_ur = -1.0;
        for (int c : frontier) {
            double cu = su + u[c], ce = se + E[c];
            double cur = (cu + ce > 0.0) ? cu / (cu + ce) : 0.0;
            if (cur > best_ur) { best_ur = cur; best = c; }
        }
        // Stop once pop_thresh is met and adding the best would drop UR below tau
        if (best == -1 || (best_ur < tau && sp >= static_cast<double>(pop_thresh))) break;
        sel_set.insert(best);
        frontier.erase(best);
        su += u[best]; se += E[best]; sp += P[best];
        for (int w : nb[best]) {
            if (!sel_set.count(w) && remaining[w]) frontier.insert(w);
        }
    }
    return std::vector<int>(sel_set.begin(), sel_set.end());
}

// -----------------------------------------------------------------------
// Full greedy-snake implementation (mirrors greedy_snake_hint)
// Phase 1: expand root first, then exhaust remaining high-UR seeds.
// Phase 2: merge groups adjacent to root that keep combined UR >= tau.
// -----------------------------------------------------------------------
static std::vector<int> greedy_snake_impl(
    const std::vector<std::vector<int>>& nb,
    const std::vector<double>& u,
    const std::vector<double>& E,
    const std::vector<double>& P,
    double tau,
    int pop_thresh,
    int root
) {
    int N = static_cast<int>(nb.size());

    // UR per node (for seed selection)
    std::vector<double> UR(N);
    for (int i = 0; i < N; ++i) {
        double d = u[i] + E[i]; UR[i] = d > 0.0 ? u[i] / d : 0.0;
    }

    std::vector<int>  assigned(N, -1);
    std::vector<bool> remaining(N, true);
    std::vector<std::vector<int>> groups;

    // Phase 1a: expand from root first
    std::vector<int> root_vec = expand_impl(root, nb, u, E, P, tau, pop_thresh, remaining);
    for (int v : root_vec) { assigned[v] = 0; remaining[v] = false; }
    groups.push_back(root_vec);

    // Phase 1b: exhaust remaining high-UR seeds
    while (true) {
        int best_seed = -1; double best_ur = -1.0;
        for (int i = 0; i < N; ++i) {
            if (remaining[i] && UR[i] >= tau && UR[i] > best_ur) {
                best_ur = UR[i]; best_seed = i;
            }
        }
        if (best_seed == -1) break;
        std::vector<int> sel = expand_impl(best_seed, nb, u, E, P, tau, pop_thresh, remaining);
        int gid = static_cast<int>(groups.size());
        for (int v : sel) { assigned[v] = gid; remaining[v] = false; }
        groups.push_back(sel);
    }

    // Phase 2: merge groups adjacent to root group that keep UR >= tau
    // Running totals for root group (updated incrementally on each merge)
    double root_u = 0.0, root_e = 0.0;
    for (int v : groups[0]) { root_u += u[v]; root_e += E[v]; }

    std::unordered_set<int> pending;
    for (int v : groups[0]) {
        for (int w : nb[v]) {
            int g = assigned[w];
            if (g > 0 && !groups[g].empty()) pending.insert(g);
        }
    }

    while (!pending.empty()) {
        int gid = *pending.begin(); pending.erase(pending.begin());
        if (groups[gid].empty()) continue;
        // Check combined UR
        double cand_u = 0.0, cand_e = 0.0;
        for (int v : groups[gid]) { cand_u += u[v]; cand_e += E[v]; }
        double combined_u = root_u + cand_u, combined_e = root_e + cand_e;
        if (combined_u + combined_e > 0.0 && combined_u / (combined_u + combined_e) >= tau) {
            // Merge
            for (int v : groups[gid]) {
                assigned[v] = 0;
                groups[0].push_back(v);
            }
            root_u = combined_u; root_e = combined_e;
            // Discover newly adjacent groups
            for (int v : groups[gid]) {
                for (int w : nb[v]) {
                    int g = assigned[w];
                    if (g > 0 && !groups[g].empty()) pending.insert(g);
                }
            }
            groups[gid].clear();
        }
    }

    std::sort(groups[0].begin(), groups[0].end());
    return groups[0];
}

// -----------------------------------------------------------------------
// Reverse-prune implementation (mirrors reverse_prune_hint)
// Start with all nodes selected; iteratively drop the node that buys the
// most UR improvement per unit of unemployment sacrificed.
// Never drops: the root, articulation points, nodes that break pop_thresh.
// -----------------------------------------------------------------------
static std::vector<int> reverse_prune_impl(
    const std::vector<std::vector<int>>& nb,
    const std::vector<double>& u,
    const std::vector<double>& E,
    const std::vector<double>& P,
    double tau,
    int pop_thresh,
    int root
) {
    int N = static_cast<int>(nb.size());
    std::vector<bool> selected(N, true);

    double U = 0.0, Ev = 0.0, Pv = 0.0;
    for (int i = 0; i < N; ++i) { U += u[i]; Ev += E[i]; Pv += P[i]; }

    while (ur_of(U, Ev) < tau) {
        std::unordered_set<int> art = art_pts_impl(nb, selected);
        double cur_ur = ur_of(U, Ev);
        int best = -1;
        double best_score = -std::numeric_limits<double>::infinity();

        for (int i = 0; i < N; ++i) {
            if (!selected[i] || i == root || art.count(i)) continue;
            if (Pv - P[i] < static_cast<double>(pop_thresh)) continue;
            double new_u = U - u[i], new_e = Ev - E[i];
            double delta_ur = ur_of(new_u, new_e) - cur_ur;
            if (delta_ur <= 0.0) continue;
            double score = (u[i] > 0.0) ? delta_ur / u[i]
                                        : std::numeric_limits<double>::infinity();
            if (score > best_score) { best_score = score; best = i; }
        }
        if (best == -1) break;
        selected[best] = false;
        U -= u[best]; Ev -= E[best]; Pv -= P[best];
    }

    std::vector<int> result;
    result.reserve(N);
    for (int i = 0; i < N; ++i) if (selected[i]) result.push_back(i);
    return result;
}

// -----------------------------------------------------------------------
// Validation helper
// -----------------------------------------------------------------------
static bool hint_valid(
    const std::vector<int>& S,
    const std::vector<double>& u,
    const std::vector<double>& E,
    const std::vector<double>& P,
    double tau,
    int pop_thresh
) {
    if (S.empty()) return false;
    double su = 0.0, se = 0.0, sp = 0.0;
    for (int i : S) { su += u[i]; se += E[i]; sp += P[i]; }
    return sp >= static_cast<double>(pop_thresh) && ur_of(su, se) >= tau;
}

static double hint_obj(const std::vector<int>& S, const std::vector<double>& u) {
    double s = 0.0; for (int i : S) s += u[i]; return s;
}

// -----------------------------------------------------------------------
// EXPORTED: cpp_articulation_points
// -----------------------------------------------------------------------
//' Find articulation points in the induced subgraph on selected nodes
//'
//' @param nb_local R List of 0-based integer neighbor vectors
//' @param selected logical vector (TRUE = node is selected / active)
//' @return 0-based integer vector of cut vertices
//' @export
// [[Rcpp::export]]
IntegerVector cpp_articulation_points(List nb_local, LogicalVector selected) {
    std::vector<std::vector<int>> nb;
    parse_nb(nb_local, nb);
    std::vector<bool> sel(selected.begin(), selected.end());
    std::unordered_set<int> art = art_pts_impl(nb, sel);
    std::vector<int> res(art.begin(), art.end());
    std::sort(res.begin(), res.end());
    return IntegerVector(res.begin(), res.end());
}

// -----------------------------------------------------------------------
// EXPORTED: cpp_greedy_snake_hint
// -----------------------------------------------------------------------
//' Greedy-snake warm-start hint for the CP-SAT ASU solver
//'
//' Mirrors \code{greedy_snake_hint()} in \file{inst/python/asu_cpsat.py}.
//' Phase 1 expands from \code{root_local} first, then exhausts remaining
//' high-UR seeds.  Phase 2 merges adjacent groups that keep UR >= tau.
//'
//' @param nb_local R List of 0-based integer neighbor vectors
//' @param u_g numeric vector of unemployed counts (same length as nb_local)
//' @param E_g numeric vector of employed counts
//' @param P_g numeric vector of population counts
//' @param tau UR threshold (e.g. 0.0645)
//' @param pop_thresh minimum population for a valid ASU
//' @param root_local 0-based index of the root/seed tract
//' @return sorted 0-based integer vector of selected tract indices
//' @export
// [[Rcpp::export]]
IntegerVector cpp_greedy_snake_hint(
    List nb_local,
    NumericVector u_g, NumericVector E_g, NumericVector P_g,
    double tau, int pop_thresh, int root_local
) {
    std::vector<std::vector<int>> nb;
    parse_nb(nb_local, nb);
    std::vector<double> u(u_g.begin(), u_g.end());
    std::vector<double> E(E_g.begin(), E_g.end());
    std::vector<double> P(P_g.begin(), P_g.end());
    std::vector<int> res = greedy_snake_impl(nb, u, E, P, tau, pop_thresh, root_local);
    return IntegerVector(res.begin(), res.end());
}

// -----------------------------------------------------------------------
// EXPORTED: cpp_reverse_prune_hint
// -----------------------------------------------------------------------
//' Reverse-prune warm-start hint for the CP-SAT ASU solver
//'
//' Mirrors \code{reverse_prune_hint()} in \file{inst/python/asu_cpsat.py}.
//' Starts with all nodes selected and iteratively drops the node that
//' buys the most UR improvement per unemployed person dropped.
//' Never drops the root, articulation points, or nodes that break pop_thresh.
//'
//' @param nb_local R List of 0-based integer neighbor vectors
//' @param u_g numeric vector of unemployed counts
//' @param E_g numeric vector of employed counts
//' @param P_g numeric vector of population counts
//' @param tau UR threshold
//' @param pop_thresh minimum population for a valid ASU
//' @param root_local 0-based index of the root/seed tract
//' @return sorted 0-based integer vector of selected tract indices
//' @export
// [[Rcpp::export]]
IntegerVector cpp_reverse_prune_hint(
    List nb_local,
    NumericVector u_g, NumericVector E_g, NumericVector P_g,
    double tau, int pop_thresh, int root_local
) {
    std::vector<std::vector<int>> nb;
    parse_nb(nb_local, nb);
    std::vector<double> u(u_g.begin(), u_g.end());
    std::vector<double> E(E_g.begin(), E_g.end());
    std::vector<double> P(P_g.begin(), P_g.end());
    std::vector<int> res = reverse_prune_impl(nb, u, E, P, tau, pop_thresh, root_local);
    return IntegerVector(res.begin(), res.end());
}

// -----------------------------------------------------------------------
// EXPORTED: cpp_compute_asu_hint
// -----------------------------------------------------------------------
//' Compute the best warm-start hint for a CP-SAT ASU window
//'
//' Runs both the greedy-snake and reverse-prune algorithms and returns the
//' one with higher total unemployment (i.e., the better objective hint).
//'
//' @param nb_local R List of 0-based integer neighbor vectors
//' @param u_g numeric vector of unemployed counts
//' @param E_g numeric vector of employed counts
//' @param P_g numeric vector of population counts
//' @param tau UR threshold
//' @param pop_thresh minimum population for a valid ASU
//' @param root_local 0-based index of the root/seed tract
//' @return Named list: \code{hint} (0-based sorted IntegerVector),
//'   \code{source} (character: "greedy_snake", "reverse_prune", or ""),
//'   \code{unemp} (integer total unemployment in hint),
//'   \code{valid} (logical: TRUE if hint meets tau and pop_thresh)
//' @export
// [[Rcpp::export]]
List cpp_compute_asu_hint(
    List nb_local,
    NumericVector u_g, NumericVector E_g, NumericVector P_g,
    double tau, int pop_thresh, int root_local
) {
    std::vector<std::vector<int>> nb;
    parse_nb(nb_local, nb);
    std::vector<double> u(u_g.begin(), u_g.end());
    std::vector<double> E(E_g.begin(), E_g.end());
    std::vector<double> P(P_g.begin(), P_g.end());

    std::vector<int> snake = greedy_snake_impl(nb, u, E, P, tau, pop_thresh, root_local);
    std::vector<int> prune = reverse_prune_impl(nb, u, E, P, tau, pop_thresh, root_local);

    bool snake_ok = hint_valid(snake, u, E, P, tau, pop_thresh);
    bool prune_ok = hint_valid(prune, u, E, P, tau, pop_thresh);

    std::vector<int> best_hint;
    std::string source;
    if (prune_ok && (!snake_ok || hint_obj(prune, u) > hint_obj(snake, u))) {
        best_hint = prune; source = "reverse_prune";
    } else if (snake_ok) {
        best_hint = snake; source = "greedy_snake";
    } else if (prune_ok) {
        best_hint = prune; source = "reverse_prune";
    } else {
        source = "";
    }

    double total_u = hint_obj(best_hint, u);
    bool valid = snake_ok || prune_ok;

    return List::create(
        Named("hint")   = IntegerVector(best_hint.begin(), best_hint.end()),
        Named("source") = source,
        Named("unemp")  = static_cast<int>(std::round(total_u)),
        Named("valid")  = valid
    );
}
