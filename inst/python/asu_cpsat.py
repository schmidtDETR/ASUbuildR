#!/usr/bin/env python3
"""
asu_cpsat.py — Build Areas of Substantial Unemployment (ASUs) with OR-Tools CP-SAT.

Modified to ensure:
1. Seeds are selected from unassigned tracts with high unemployment rate
2. Stops when no remaining unassigned tract has UR >= tau (6.45%)

Two ways to provide adjacency (contiguity):
  A) Precomputed neighbors JSON (recommended on servers without a Geo stack)
     - JSON format: list of lists of tract indices (0- or 1-based). E.g., [[1,4],[2],...]
  B) Geometry file (GeoPackage/Shapefile) + libpysal Queen contiguity

Input table must contain:
  geoid, tract_ASU_unemp, tract_ASU_emp, tract_pop2024

Examples:
  # Using precomputed neighbors (fastest to deploy)
  python asu_cpsat.py --input OR_asu26.xlsx --neighbors nb_queen.json \
    --tau 0.0645 --pop-thresh 10000 --max-asus 30 --time-limit 1200 --workers 8 --verbose

  # Compute queen contiguity from geometry (requires geopandas, libpysal)
  python asu_cpsat.py --input OR_asu26.xlsx --geometry tracts_or_2024.gpkg \
    --geom-col geometry --geoid-col GEOID \
    --tau 0.0645 --pop-thresh 10000 --max-asus 30 --time-limit 1200 --workers 8 --verbose
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import time
from typing import List, Optional, Dict, Tuple

import numpy as np
import pandas as pd
from ortools.sat.python import cp_model

# Optional (only needed if you compute contiguity on the fly)
try:
    import geopandas as gpd
    from libpysal.weights import Queen
    from shapely.validation import make_valid as shapely_make_valid
except Exception:
    gpd = None
    Queen = None
    shapely_make_valid = None


# ---------- Helpers ----------
def as_fraction_tau(tau: float) -> Tuple[int, int]:
    """Represent k = tau/(1-tau) as num/den using exact integers when tau has 4 decimals."""
    T = int(round(tau * 10000))
    one_minus = 10000 - T
    g = math.gcd(T, one_minus)
    return T // g, one_minus // g  # k = num/den


def ur_of(u_sum: int, E_sum: int) -> float:
    return 0.0 if (u_sum + E_sum) == 0 else u_sum / (u_sum + E_sum)


def bfs_ball(nb: List[List[int]], center: int, r: int, allowed: np.ndarray) -> List[int]:
    allowed_set = set(int(a) for a in allowed)
    vis = {center}
    frontier = [center]
    layer = 0
    while layer < r:
        nxt = []
        for v in frontier:
            for w in nb[v]:
                if (w in allowed_set) and (w not in vis):
                    vis.add(w)
                    nxt.append(w)
        if not nxt:
            break
        frontier = list(set(nxt))
        layer += 1
    return sorted(vis)


def greedy_snake_hint(
    nb_local: List[List[int]],
    u_g: np.ndarray, # unemployment counts per tract
    E_g: np.ndarray, # employment counts per tract
    P_g: np.ndarray, # population counts per tract
    tau: float, # unemployment rate threshold for ASU selection
    pop_thresh: int, # minimum population threshold for ASU selection
    root_local: int,
) -> List[int]:
    """
    Run Simple Snake across the full local window, then combine touching groups.
    Phase 1: seed root_local first, then exhaust all remaining high-UR seeds.
    Phase 2: merge any groups adjacent to root's group that keep UR >= tau.
    Returns the resulting merged group (local indices) as the CP-SAT warm-start hint.
    """
    N = len(nb_local)
    UR = u_g / np.maximum(u_g + E_g, 1e-12)
    assigned = np.full(N, -1, dtype=int)   # group id per node, -1 = unassigned
    remaining = np.ones(N, dtype=bool)

    def _expand(seed: int) -> set:
        sel = {seed}
        sel_u = int(u_g[seed])
        sel_e = int(E_g[seed])
        sel_p = int(P_g[seed])
        # incremental frontier: O(N·degree) total instead of O(N²·degree)
        frontier = {w for w in nb_local[seed] if remaining[w]}
        while frontier:
            best, best_ur = None, -1.0
            for cand in frontier:
                cu = sel_u + int(u_g[cand])
                ce = sel_e + int(E_g[cand])
                cur = cu / (cu + ce) if cu + ce > 0 else 0.0
                if cur > best_ur:
                    best_ur, best = cur, cand
            if best is None or (best_ur < tau and sel_p >= pop_thresh):
                break
            sel.add(best)
            frontier.discard(best)
            sel_u += int(u_g[best])
            sel_e += int(E_g[best])
            sel_p += int(P_g[best])
            for w in nb_local[best]:
                if w not in sel and remaining[w]:
                    frontier.add(w)
        return sel

    # Phase 1: seed root_local first, then exhaust remaining high-UR seeds
    groups: List[set] = []

    root_set = _expand(root_local)
    root_gid = 0
    for v in root_set:
        assigned[v] = root_gid
        remaining[v] = False
    groups.append(root_set)

    while True:
        rem = np.where(remaining & (UR >= tau))[0]
        if rem.size == 0:
            break
        seed = int(rem[np.argmax(UR[rem])])
        sel = _expand(seed)
        gid = len(groups)
        for v in sel:
            assigned[v] = gid
            remaining[v] = False
        groups.append(sel)

    # Phase 2: merge groups that touch root_gid and keep combined UR >= tau
    pending: set = set()
    for v in groups[root_gid]:
        for w in nb_local[v]:
            g = assigned[w]
            if g not in (-1, root_gid) and len(groups[g]) > 0:
                pending.add(g)

    while pending:
        gid = pending.pop()
        if len(groups[gid]) == 0:
            continue
        combined = list(groups[root_gid] | groups[gid])
        cu = int(u_g[combined].sum())
        ce = int(E_g[combined].sum())
        if cu + ce > 0 and cu / (cu + ce) >= tau:
            new_members = groups[gid]
            for v in new_members:
                assigned[v] = root_gid
            groups[root_gid] = groups[root_gid] | new_members
            groups[gid] = set()
            # discover neighbours newly exposed by the merged members
            for v in new_members:
                for w in nb_local[v]:
                    g = assigned[w]
                    if g not in (-1, root_gid) and len(groups[g]) > 0:
                        pending.add(g)

    return sorted(groups[root_gid])


def _articulation_points(nb_local: List[List[int]], selected: np.ndarray) -> set:
    """
    Iterative Tarjan articulation-point finder restricted to the induced subgraph
    on `selected` nodes. A cut vertex's removal disconnects the remainder of its
    connected component, so these are never valid drop candidates for
    reverse_prune_hint. Mirrors the iterative low-link style of
    _bridge_edge_bounds but for vertices instead of edges.
    Returns a set of local node indices that are cut vertices.
    """
    N = len(nb_local)
    disc = [-1] * N
    low = [0] * N
    parent = [-1] * N
    skipped_parent = [False] * N
    root_children = [0] * N
    is_art: set = set()
    timer = 0

    for start in range(N):
        if not selected[start] or disc[start] != -1:
            continue
        stack = [(start, iter(nb_local[start]))]
        disc[start] = low[start] = timer
        timer += 1
        while stack:
            u, it = stack[-1]
            recursed = False
            for w in it:
                if not selected[w]:
                    continue
                if w == parent[u] and not skipped_parent[u]:
                    skipped_parent[u] = True
                    continue
                if disc[w] == -1:
                    parent[w] = u
                    if parent[u] == -1:
                        root_children[start] += 1
                    disc[w] = low[w] = timer
                    timer += 1
                    stack.append((w, iter(nb_local[w])))
                    recursed = True
                    break
                else:
                    low[u] = min(low[u], disc[w])
            if not recursed:
                stack.pop()
                if stack:
                    p = stack[-1][0]
                    low[p] = min(low[p], low[u])
                    if parent[p] != -1 and low[u] >= disc[p]:
                        is_art.add(p)
        if root_children[start] > 1:
            is_art.add(start)
    return is_art


def reverse_prune_hint(
    nb_local: List[List[int]],
    u_g: np.ndarray,
    E_g: np.ndarray,
    P_g: np.ndarray,
    tau: float,
    pop_thresh: int,
    root_local: int,
) -> List[int]:
    """
    Warm start via reverse pruning: start with every tract in the window selected
    (nb_local is itself a connected BFS ball, so this is guaranteed connected),
    then repeatedly drop the tract that buys the most UR improvement per unit of
    unemployment sacrificed, scored as delta_ur / u_dropped (tracts that raise UR
    at zero unemployment cost score as +inf and are dropped first). Never drops
    the root or a current cut vertex of the selected induced subgraph (that would
    disconnect the remainder), and never drops below pop_thresh. Stops once
    UR >= tau (success) or no valid drop remains (failure -- caller validates via
    component_ok and falls back to other warm-start sources if so).
    """
    N = len(nb_local)
    selected = np.ones(N, dtype=bool)
    U_sum = int(u_g.sum())
    E_sum = int(E_g.sum())
    P_sum = int(P_g.sum())

    while ur_of(U_sum, E_sum) < tau:
        cur_ur = ur_of(U_sum, E_sum)
        cut_vertices = _articulation_points(nb_local, selected)

        new_u = U_sum - u_g
        new_e = E_sum - E_g
        new_ur = new_u / np.maximum(new_u + new_e, 1e-12)
        delta_ur = new_ur - cur_ur

        droppable = selected.copy()
        droppable[root_local] = False
        for v in cut_vertices:
            droppable[v] = False
        droppable &= delta_ur > 0
        droppable &= (P_sum - P_g) >= pop_thresh

        cand_idx = np.where(droppable)[0]
        if cand_idx.size == 0:
            break  # can't reach tau this way without disconnecting or breaking pop_thresh

        u_drop = u_g[cand_idx]
        score = np.where(u_drop > 0, delta_ur[cand_idx] / np.maximum(u_drop, 1), np.inf)
        best = int(cand_idx[np.argmax(score)])

        selected[best] = False
        U_sum -= int(u_g[best])
        E_sum -= int(E_g[best])
        P_sum -= int(P_g[best])

    return sorted(int(v) for v in range(N) if selected[v])


def _spanning_tree_flows(hint: List[int], nb_local: List[List[int]], root_local: int) -> Dict[Tuple[int, int], int]:
    """
    Compute single-commodity flow values for the BFS spanning tree of hint.
    Each tree edge (parent→child) carries flow = subtree size at child.
    Returns a dict of (i, j) -> flow for non-zero entries only.
    """
    if not hint:
        return {}
    hint_set = set(hint)
    parent: Dict[int, Optional[int]] = {root_local: None}
    children: Dict[int, List[int]] = {v: [] for v in hint}
    queue = [root_local]
    order: List[int] = []
    while queue:
        v = queue.pop(0)
        order.append(v)
        for w in nb_local[v]:
            if w in hint_set and w not in parent:
                parent[w] = v
                children[v].append(w)
                queue.append(w)
    size: Dict[int, int] = {v: 1 for v in hint}
    for v in reversed(order):
        for c in children[v]:
            size[v] += size[c]
    flows: Dict[Tuple[int, int], int] = {}
    for v in hint:
        p = parent[v]
        if p is not None:
            flows[(p, v)] = size[v]
    return flows


def component_ok(S: List[int], u: np.ndarray, E: np.ndarray, P: np.ndarray,
                 tau: float, pop_thresh: int, nb: List[List[int]]) -> bool:
    if not S:
        return False
    Sset = set(S)
    # connectivity (BFS)
    seen = {S[0]}
    Q = [S[0]]
    while Q:
        v = Q.pop()
        for w in nb[v]:
            if (w in Sset) and (w not in seen):
                seen.add(w)
                Q.append(w)
    if len(seen) != len(S):
        return False
    su, sE, sP = int(u[S].sum()), int(E[S].sum()), int(P[S].sum())
    return (sP >= pop_thresh) and (ur_of(su, sE) >= tau)


def can_hit_tau(u: np.ndarray, E: np.ndarray, P: np.ndarray,
                nb_local: List[List[int]], tau: float, pop_thresh: int) -> bool:
    """Quick optimistic screen: if no component can meet UR/pop, skip solving."""
    if len(u) == 0:
        return False
    UR = u / np.maximum(u + E, 1e-12)
    if UR.max(initial=0.0) < tau:
        return False
    # BFS-ball windows are always connected; single-component knapsack is sufficient
    num, den = as_fraction_tau(tau)
    D = den * u - num * E
    rho = D / np.maximum(P, 1e-12)
    ord_idx = np.argsort(-rho)
    need = pop_thresh
    cumD = 0.0
    for j in ord_idx:
        pj, dj = int(P[j]), float(D[j])
        if pj <= 0:
            continue
        take = min(pj, need)
        cumD += dj * (take / pj)
        need -= take
        if need <= 0:
            break
    return need <= 0 and cumD >= 0


def queen_neighbors_from_geometries(gdf: "gpd.GeoDataFrame", geom_col: str = "geometry") -> List[List[int]]:
    if gpd is None or Queen is None:
        raise RuntimeError("geopandas + libpysal required to compute contiguity from geometry.")
    if geom_col not in gdf.columns:
        raise ValueError(f"Geometry column '{geom_col}' not found.")

    gdf = gdf.reset_index(drop=True)
    # basic validity repair
    if hasattr(gdf.geometry, "is_valid"):
        invalid = ~gdf.geometry.is_valid
        if invalid.any():
            gdf.loc[invalid, geom_col] = gdf.loc[invalid, geom_col].buffer(0)
            invalid = ~gdf.geometry.is_valid
            if invalid.any() and shapely_make_valid is not None:
                gdf.loc[invalid, geom_col] = gdf.loc[invalid, geom_col].apply(shapely_make_valid)

    W = Queen.from_dataframe(gdf, ids=list(range(len(gdf))))
    nb = [[] for _ in range(len(gdf))]
    for i, neigh in W.neighbors.items():
        nb[i] = sorted(neigh)
    return nb


def contract_high_ur_nodes(
    nb_local: List[List[int]],
    u_g: np.ndarray,
    E_g: np.ndarray,
    P_g: np.ndarray,
    tau: float,
) -> Tuple[List[List[int]], np.ndarray, np.ndarray, np.ndarray, List[List[int]], np.ndarray]:
    """
    Fuse each connected cluster of UR>=tau tracts into one super-node (exact, lossless).
    Mediant inequality: mixing two UR>=tau distributions keeps combined UR>=tau, so
    every cluster member co-occurs in every optimal solution — no branching needed.
    Returns: nb_r, u_r, E_r, P_r, expand, node_map
      expand[ri]  = sorted list of original local indices for reduced node ri
      node_map[v] = reduced index for original local node v
    """
    N = len(nb_local)
    num, den = as_fraction_tau(tau)
    high_set = {
        i for i in range(N)
        if int(den) * int(u_g[i]) - int(num) * int(E_g[i]) >= 0
    }

    # BFS over induced high-UR subgraph to label connected components
    comp = np.full(N, -1, dtype=int)
    c = 0
    for v in sorted(high_set):
        if comp[v] >= 0:
            continue
        stk = [v]
        comp[v] = c
        while stk:
            cur = stk.pop()
            for w in nb_local[cur]:
                if w in high_set and comp[w] < 0:
                    comp[w] = c
                    stk.append(w)
        c += 1
    n_super = c

    # Reduced indices: super-nodes 0..n_super-1, then one index per low-UR node
    node_map = np.empty(N, dtype=int)
    next_low = n_super
    for v in range(N):
        if v in high_set:
            node_map[v] = comp[v]
        else:
            node_map[v] = next_low
            next_low += 1
    N_r = next_low

    expand: List[List[int]] = [[] for _ in range(N_r)]
    for v in range(N):
        expand[node_map[v]].append(v)

    u_r = np.array([int(u_g[expand[ri]].sum()) for ri in range(N_r)], dtype=np.int64)
    E_r = np.array([int(E_g[expand[ri]].sum()) for ri in range(N_r)], dtype=np.int64)
    P_r = np.array([int(P_g[expand[ri]].sum()) for ri in range(N_r)], dtype=np.int64)

    adj_sets: List[set] = [set() for _ in range(N_r)]
    for v in range(N):
        rv = int(node_map[v])
        for w in nb_local[v]:
            rw = int(node_map[w])
            if rv != rw:
                adj_sets[rv].add(rw)
                adj_sets[rw].add(rv)
    nb_r = [sorted(s) for s in adj_sets]

    return nb_r, u_r, E_r, P_r, expand, node_map


# ---------- CP-SAT core: solve one ASU within a window ----------
def _bridge_edge_bounds(nb_local: List[List[int]]) -> Dict[Tuple[int, int], int]:
    """
    Find bridges (cut edges) of the FIXED underlying graph and, for each, the size
    of each side's component when that bridge is removed. Sound (unlike a
    reference-spanning-tree bound of the *selected* subgraph, which depends on
    the unknown selection and was proven unsound this session): a bridge splits
    the known, fixed graph into exactly two components regardless of what ends up
    selected, so the flow crossing it in either direction can never exceed that
    side's node count. Non-bridge edges (on any cycle) get no tightening here --
    flow could route around the cycle, so no valid per-edge bound below N-1 can
    be derived from local structure alone.
    Returns {(u, v): bound} for every directed bridge edge in both directions.
    """
    N = len(nb_local)
    disc = [-1] * N
    low = [0] * N
    subtree_size = [1] * N
    parent = [-1] * N
    skipped_parent = [False] * N
    bounds: Dict[Tuple[int, int], int] = {}
    timer = 0

    for start in range(N):
        if disc[start] != -1:
            continue
        comp_nodes = 0
        comp_bridges: List[Tuple[int, int]] = []
        stack = [(start, iter(nb_local[start]))]
        disc[start] = low[start] = timer
        timer += 1
        while stack:
            u, it = stack[-1]
            recursed = False
            for w in it:
                if w == parent[u] and not skipped_parent[u]:
                    skipped_parent[u] = True
                    continue
                if disc[w] == -1:
                    parent[w] = u
                    disc[w] = low[w] = timer
                    timer += 1
                    stack.append((w, iter(nb_local[w])))
                    recursed = True
                    break
                else:
                    low[u] = min(low[u], disc[w])
            if not recursed:
                stack.pop()
                comp_nodes += 1
                if stack:
                    p = stack[-1][0]
                    low[p] = min(low[p], low[u])
                    subtree_size[p] += subtree_size[u]
                    if low[u] > disc[p]:
                        # (p, u) is a bridge -- removing it isolates u's DFS
                        # subtree from the rest of this connected component.
                        # subtree_size[u] is already final (u's own subtree is
                        # fully explored before u is popped), but the *other*
                        # side's size needs the component's final total node
                        # count, which isn't known until this whole DFS
                        # finishes -- so defer and fill in below.
                        bounds[(p, u)] = subtree_size[u]
                        comp_bridges.append((u, p))
        for (u, p) in comp_bridges:
            bounds[(u, p)] = comp_nodes - subtree_size[u]
    return bounds


class CpsatResult:
    def __init__(self, sel_idx_local: List[int], root_local: int, obj: int, status: str):
        self.sel_idx_local = sel_idx_local
        self.root_local = root_local
        self.obj = obj
        self.status = status


def solve_one_asu_cpsat(
    nb_local: List[List[int]],
    u_g: np.ndarray,
    E_g: np.ndarray,
    P_g: np.ndarray,
    tau: float,
    pop_thresh: int,
    root_local: int,
    time_limit: int = 1200,
    workers: int = 8,
    log: bool = True,
    rel_gap: Optional[float] = None,
    hint: Optional[List[int]] = None,
    hint_obj: Optional[int] = None,
    forced_selected: Optional[List[int]] = None,
    cluster_groups: Optional[List[List[int]]] = None,
) -> Optional[CpsatResult]:
    """
        Connectivity via iterative vertex-separator cuts. Each disconnected incumbent
        adds valid constraints requiring a selected path from its components to the root.
    Objective: maximize Σ u_i x_i
    """
    N = len(nb_local)
    if N == 0:
        return None

    model = cp_model.CpModel()

    # Decision variables
    x = [model.NewBoolVar(f"x_{i}") for i in range(N)]

    # Selecting the root permits its entire connected high-UR component at no cost.
    forced_set = set(forced_selected or [root_local])
    forced_set.add(root_local)
    for i in forced_set:
        model.Add(x[i] == 1)

    # Valid for ANY connected subgraph: every selected non-root node must have at
    # least one selected neighbor (it can't be reached from root otherwise). Free
    # to add upfront (no separation needed) and tightens both the cut phase and the
    # flow phase since they share this same model.
    for v in range(N):
        if v == root_local:
            continue
        if nb_local[v]:
            model.Add(x[v] <= sum(x[w] for w in nb_local[v]))
        else:
            model.Add(x[v] == 0)

    # Tie every connected UR>=tau cluster's members together (x[v] == x[anchor]).
    # Sound for a maximization objective: since each member's own UR ratio is >=
    # tau, mixing in an adjacent cluster member via the mediant inequality can only
    # keep the combined ratio >= tau, population only grows, and connectivity is
    # preserved (cluster is internally connected) -- so once any member is
    # selected, including the rest is always weakly better. This reduces the
    # solver's effective search space without touching the flow graph at all (the
    # graph-based contraction that previously caused presolve infeasibility is
    # NOT used here -- nb_local/edges stay exactly as-is).
    if cluster_groups:
        for group in cluster_groups:
            if len(group) > 1:
                anchor = group[0]
                for v in group[1:]:
                    model.Add(x[v] == x[anchor])

    # Population threshold
    pop_expr = sum(int(P_g[i]) * x[i] for i in range(N))
    model.Add(pop_expr >= int(pop_thresh))

    # UR >= tau as exact integer linear inequality
    num, den = as_fraction_tau(tau)
    lhs = sum(int(den) * int(u_g[i]) * x[i] for i in range(N)) \
        - sum(int(num) * int(E_g[i]) * x[i] for i in range(N))
    model.Add(lhs >= 0)

    # Objective: maximize unemployment captured
    obj_expr = sum(int(u_g[i]) * x[i] for i in range(N))
    model.Maximize(obj_expr)

    # Warm-start with the connected greedy snake solution.
    if hint is not None:
        hint_set = set(hint)
        for i in range(N):
            model.AddHint(x[i], 1 if i in hint_set else 0)

    # Lower bound: reject solutions worse than the greedy warm start
    if hint_obj is not None and hint_obj > 0:
        model.Add(obj_expr >= hint_obj)

    # Solver params
    best_connected = sorted(hint) if hint and component_ok(
        hint, u_g, E_g, P_g, tau, pop_thresh, nb_local
    ) else None
    best_obj = int(u_g[best_connected].sum()) if best_connected else -1
    lower_bound = hint_obj if (hint_obj is not None and hint_obj > 0) else -1
    start_time = time.monotonic()
    cut_round = 0
    # NOTE: a tight objective/bound in the cut-only (disconnected) relaxation does
    # NOT imply cuts are close to finding a *connected* solution -- the "price of
    # connectivity" gap can be large and take many more rounds to close than the
    # relaxation's bound suggests (measured: 15 rounds over ~30s shrank components
    # from 18->8 without ever reaching a single connected component). The exact
    # flow-based phase is what actually *guarantees* progress toward a connected
    # answer, so it must keep the majority of the time budget; cuts are only a
    # cheap pre-pass to prune obviously-disconnected structure.
    #
    # NOTE: tried removing this cut pre-pass entirely (going straight to the
    # exact flow model) and A/B tested it on real Colorado data -- it was a
    # clear regression (0.35%->0.79% gap, 75,357->75,204 unemp @300s) with no
    # wall-time savings (still used the full 300s budget either way). The
    # boundary constraints these cuts add evidently still prune the flow
    # phase's search space usefully even when they never converge to a single
    # connected component. See SKILL.md.
    cut_time_budget = min(20.0, max(2.0, float(time_limit) * 0.15))
    stall_rounds = 0
    prev_num_components: Optional[int] = None
    first_components: Optional[int] = None

    while cut_round < 15:
        elapsed = time.monotonic() - start_time
        remaining_for_cuts = cut_time_budget - elapsed
        if remaining_for_cuts <= 0:
            break

        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = max(1, int(workers))
        solver.parameters.max_time_in_seconds = remaining_for_cuts
        solver.parameters.log_search_progress = False  # silent; summary logged after loop
        solver.parameters.cp_model_presolve = True
        solver.parameters.linearization_level = 2

        status = solver.Solve(model)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            break

        selected = [i for i in range(N) if solver.BooleanValue(x[i])]
        selected_set = set(selected)
        root_component = {root_local}
        stack = [root_local]
        while stack:
            v = stack.pop()
            for w in nb_local[v]:
                if w in selected_set and w not in root_component:
                    root_component.add(w)
                    stack.append(w)

        if len(root_component) == len(selected):
            objective = int(round(solver.ObjectiveValue()))
            if objective > best_obj:
                best_connected, best_obj = selected, objective
                if best_obj > lower_bound:
                    model.Add(obj_expr >= best_obj)
                    lower_bound = best_obj
            if status == cp_model.OPTIMAL:
                return CpsatResult(selected, root_local, objective, "OPTIMAL")
            break

        unseen = selected_set - root_component
        components: List[set] = []
        while unseen:
            seed = unseen.pop()
            component = {seed}
            stack = [seed]
            while stack:
                v = stack.pop()
                for w in nb_local[v]:
                    if w in unseen:
                        unseen.remove(w)
                        component.add(w)
                        stack.append(w)
            components.append(component)

        if first_components is None:
            first_components = len(components)
        if prev_num_components is not None and len(components) >= prev_num_components:
            stall_rounds += 1
        else:
            stall_rounds = 0
        prev_num_components = len(components)

        for component in components:
            boundary = {
                w for v in component for w in nb_local[v]
                if w not in component
            }
            if boundary:
                boundary_expr = sum(x[w] for w in boundary)
                for v in component:
                    model.Add(x[v] <= boundary_expr)
            else:
                for v in component:
                    model.Add(x[v] == 0)
        cut_round += 1
        if stall_rounds >= 3:
            break

    # Finish with exact single-commodity flow connectivity, strengthened by the cuts.
    edges = list(dict.fromkeys(
        (i, j) for i, neighbors in enumerate(nb_local) for j in neighbors if i != j
    ))

    # NOTE: a per-edge bound derived from a single fixed reference spanning tree
    # (e.g. subtree size) is UNSOUND on graphs with cycles -- the actual flow can
    # legitimately need to route around a different topology than any one fixed
    # tree, and a fixed-tree bound can wrongly reject genuinely feasible connected
    # selections. Verified empirically: a 5-cycle counterexample where excluding
    # one low-value node forces routing 3 units through what a BFS tree treats as
    # a capacity-2 edge. The uniform bound below is the correct, universally valid
    # one (flow on any edge can never exceed total selected nodes - 1).
    M = max(1, N - 1)
    # NOTE: _bridge_edge_bounds() gives a provably sound tighter per-edge cap
    # (validated against 400 brute-force instances + explicit counterexamples)
    # but was A/B tested on real Colorado data and was a clear regression
    # (0.52%->2.62% gap, 75,214->74,160 unemp @300s) -- see SKILL.md. Kept
    # available but NOT wired in by default.
    edge_bounds = [M] * len(edges)

    f = [model.NewIntVar(0, edge_bounds[idx], f"f_{i}_{j}") for idx, (i, j) in enumerate(edges)]
    in_edges_for = [[] for _ in range(N)]
    out_edges_for = [[] for _ in range(N)]
    for edge_index, (i, j) in enumerate(edges):
        out_edges_for[i].append(edge_index)
        in_edges_for[j].append(edge_index)
        model.Add(f[edge_index] <= edge_bounds[edge_index] * x[i])
        model.Add(f[edge_index] <= edge_bounds[edge_index] * x[j])

    selected_count = sum(x)
    for i in range(N):
        inflow = sum(f[e] for e in in_edges_for[i]) if in_edges_for[i] else 0
        outflow = sum(f[e] for e in out_edges_for[i]) if out_edges_for[i] else 0
        if i == root_local:
            model.Add(outflow - inflow == selected_count - 1)
        else:
            model.Add(inflow - outflow == x[i])

    # Warm-start flows from the best connected incumbent found so far (the cut
    # phase may have improved on the caller-supplied hint).
    flow_source = best_connected if best_connected is not None else hint
    if flow_source is not None:
        flow_hints = _spanning_tree_flows(flow_source, nb_local, root_local)
        for edge_index, edge in enumerate(edges):
            model.AddHint(f[edge_index], flow_hints.get(edge, 0))

    remaining_time = float(time_limit) - (time.monotonic() - start_time)
    if log and cut_round > 0:
        _fc = first_components or "?"
        _lc = prev_num_components if prev_num_components is not None else "?"
        print(f"  cut phase: {cut_round} round(s), components {_fc}->{_lc}, "
              f"{time_limit - remaining_time:.1f}s used; {remaining_time:.1f}s for flow phase", flush=True)
    if remaining_time <= 0:
        return CpsatResult(best_connected, root_local, best_obj, "FEASIBLE") if best_connected else None

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = max(1, int(workers))
    solver.parameters.max_time_in_seconds = remaining_time
    solver.parameters.log_search_progress = bool(log)
    solver.parameters.cp_model_presolve = True
    solver.parameters.linearization_level = 2
    # pseudo_costs/reduced_costs/default_lp/quick_restart never produced a
    # winning incumbent in verbose benchmark logs and are safe to drop.
    # NOTE: do NOT also glob-exclude "probing*"/"objective_lb_search*" --
    # A/B tested 2026-07-30 and it was a clear regression (0.52%->1.86% gap,
    # 75,214->74,623 unemp @300s) despite those workers never winning
    # directly themselves, likely via shared clauses/bounds helping LNS/LS.
    solver.parameters.ignore_subsolvers.extend([
        "pseudo_costs", "reduced_costs", "default_lp", "quick_restart",
    ])
    if rel_gap is not None:
        solver.parameters.relative_gap_limit = float(rel_gap)
    status = solver.Solve(model)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        selected = [i for i in range(N) if solver.BooleanValue(x[i])]
        objective = int(round(solver.ObjectiveValue()))
        return CpsatResult(selected, root_local, objective, solver.StatusName(status))
    if best_connected is not None:
        return CpsatResult(best_connected, root_local, best_obj, "FEASIBLE")
    return None


# ---------- Simple improver (local trades) ----------
def frontier_candidates(S: List[int], nb: List[List[int]], allowed: np.ndarray) -> List[int]:
    Sset = set(S)
    allowed_set = set(int(a) for a in allowed)
    cand = set()
    for v in S:
        for w in nb[v]:
            if (w not in Sset) and (w in allowed_set):
                cand.add(w)
    return sorted(cand)


def improve_by_trades(S0: List[int], u: np.ndarray, E: np.ndarray, P: np.ndarray, nb: List[List[int]],
                      tau: float, pop_thresh: int, allowed: np.ndarray, max_iter: int = 200) -> List[int]:
    S = sorted(set(S0))
    selected = set(S)
    sum_u = int(u[S].sum())
    sum_E = int(E[S].sum())
    sum_P = int(P[S].sum())
    for _ in range(max_iter):
        improved = False
        # Greedy adds: frontier nodes are adjacent to S so S∪{t} is always connected
        for t in sorted(frontier_candidates(S, nb, allowed), key=lambda i: u[i], reverse=True):
            next_u = sum_u + int(u[t])
            next_E = sum_E + int(E[t])
            next_P = sum_P + int(P[t])
            if next_P >= pop_thresh and ur_of(next_u, next_E) >= tau:
                selected.add(t)
                S = sorted(selected)
                sum_u, sum_E, sum_P = next_u, next_E, next_P
                improved = True
                break
        if improved:
            continue
        # Swap: drop worst u, add best neighbor
        if len(S) > 1:
            for r in sorted(S, key=lambda i: u[i]):
                reduced_u = sum_u - int(u[r])
                reduced_E = sum_E - int(E[r])
                reduced_P = sum_P - int(P[r])
                if reduced_P < pop_thresh or ur_of(reduced_u, reduced_E) < tau:
                    continue
                S2 = sorted(selected - {r})
                if not component_ok(S2, u, E, P, tau, pop_thresh, nb):
                    continue
                for a in sorted(frontier_candidates(S2, nb, allowed), key=lambda i: u[i], reverse=True):
                    next_u = reduced_u + int(u[a])
                    next_E = reduced_E + int(E[a])
                    next_P = reduced_P + int(P[a])
                    if (next_u > sum_u and next_P >= pop_thresh
                            and ur_of(next_u, next_E) >= tau):
                        selected = set(S2)
                        selected.add(a)
                        S = sorted(selected)
                        sum_u, sum_E, sum_P = next_u, next_E, next_P
                        improved = True
                        break
                if improved:
                    break
        if not improved:
            break
    return S


def _prepare_window_hint(
    nb_local: List[List[int]], u_g: np.ndarray, E_g: np.ndarray, P_g: np.ndarray,
    tau: float, pop_thresh: int, root_local: int,
) -> Dict:
    """
    Contract UR>=tau clusters, then build two candidate warm starts on the
    reduced graph -- a forward-growing greedy snake and a reverse prune that
    starts from the whole window and drops tracts down to tau -- refine each on
    the original graph, and keep whichever reaches the higher unemployment
    objective.
    """
    nb_r, u_r, E_r, P_r, expand_r, node_map_r = contract_high_ur_nodes(nb_local, u_g, E_g, P_g, tau)
    root_r = int(node_map_r[root_local])
    root_component = expand_r[root_r]
    all_local = np.arange(len(nb_local))

    def _refine(hint_r: List[int]) -> Dict:
        hint_expanded = sorted({orig for ri in hint_r for orig in expand_r[ri]})
        hint_improved = improve_by_trades(hint_expanded, u_g, E_g, P_g, nb_local, tau, pop_thresh, all_local, max_iter=100)
        hint_valid = component_ok(hint_improved, u_g, E_g, P_g, tau, pop_thresh, nb_local)
        hint_obj_val = int(u_g[hint_improved].sum()) if hint_valid else None
        return {"hint_improved": hint_improved, "hint_valid": hint_valid, "hint_obj_val": hint_obj_val}

    snake = _refine(greedy_snake_hint(nb_r, u_r, E_r, P_r, tau, pop_thresh, root_r))
    prune = _refine(reverse_prune_hint(nb_r, u_r, E_r, P_r, tau, pop_thresh, root_r))

    if prune["hint_valid"] and (not snake["hint_valid"] or prune["hint_obj_val"] > snake["hint_obj_val"]):
        best, hint_source = prune, "reverse_prune"
    else:
        best, hint_source = snake, "greedy_snake"

    return {
        "root_component": root_component,
        "n_contracted": len(nb_r),
        "hint_improved": best["hint_improved"],
        "hint_valid": best["hint_valid"],
        "hint_obj_val": best["hint_obj_val"],
        "hint_source": hint_source,
        "cluster_groups": [group for group in expand_r if len(group) > 1],
    }


# ---------- High-level multi-ASU builder ----------
def _export_window_comparison(
    w: Dict,
    sol_local: Optional[List[int]],
    df: pd.DataFrame,
    export_dir: str,
    asu_num: int,
) -> None:
    """Write tract-comparison + neighbor-list Excel for one solved window."""
    import os
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        print("  [EXPORT] skipped: openpyxl not installed (pip install openpyxl)", flush=True)
        return

    sub = w["sub"]
    nb_local: List[List[int]] = w["nb_local"]
    N = len(sub)
    geoids = [str(df.iloc[int(sub[i])]["geoid"]) for i in range(N)]
    unemp = w["u_g"]
    emp   = w["E_g"]
    pop   = w["P_g"]

    hint_set = set(w.get("hint_improved") or [])
    sol_set  = set(sol_local or [])
    root     = w["root_local"]
    hsrc     = w.get("hint_source", "")

    tract_rows = []
    for i in range(N):
        u_i, e_i, p_i = int(unemp[i]), int(emp[i]), int(pop[i])
        tract_rows.append({
            "local_idx":   i,
            "geoid":       geoids[i],
            "unemp":       u_i,
            "emp":         e_i,
            "pop":         p_i,
            "ur_pct":      round(u_i / max(u_i + e_i, 1) * 100, 4),
            "in_hint":     i in hint_set,
            "hint_source": hsrc if i in hint_set else "",
            "in_solution": i in sol_set,
            "is_root":     i == root,
        })

    seen: set = set()
    edge_rows = []
    for i in range(N):
        for j in nb_local[i]:
            edge = (min(i, j), max(i, j))
            if edge not in seen:
                seen.add(edge)
                edge_rows.append({
                    "from_idx":   i,
                    "to_idx":     j,
                    "from_geoid": geoids[i],
                    "to_geoid":   geoids[j],
                })

    os.makedirs(export_dir, exist_ok=True)
    path = os.path.join(export_dir, f"asu_{asu_num:03d}_seed{w['seed']}.xlsx")
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(tract_rows).to_excel(writer, sheet_name="tracts",    index=False)
        pd.DataFrame(edge_rows).to_excel( writer, sheet_name="neighbors", index=False)
    print(f"  [EXPORT] {path}", flush=True)


def build_many_asus_cpsat(
    df: pd.DataFrame,
    nb: List[List[int]],
    tau: float,
    pop_thresh: int,
    max_asus: int = 25,
    r_start: int = 50,
    r_step: int = 1,
    r_max: int = 50,
    hard_cap_nodes: int = 10000,
    min_pop_margin: float = 1.0,
    time_limit: int = 1200,
    workers: int = 8,
    rel_gap: Optional[float] = None,
    verbose: bool = True,
    parallel_asus: int = 4,
    merge_adjacent: bool = True,
    export_dir: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """
    Build ASUs in batches of up to `parallel_asus` disjoint candidate windows, solved
    concurrently. Two ASUs built in the same batch that end up touching (share a
    queen-contiguity edge) are merged into one: the mediant inequality guarantees
    that combining two groups whose UR is each >= tau keeps the combined UR >= tau
    (the combined ratio is a weighted average of the two, so it can't fall below
    the smaller one), and population/connectivity only improve on union.
    """
    u = df["tract_ASU_unemp"].to_numpy(dtype=np.int64)
    E = df["tract_ASU_emp"].to_numpy(dtype=np.int64)
    P = df["tract_pop2024"].to_numpy(dtype=np.int64)
    UR = u / np.maximum(u + E, 1e-12)

    n = len(df)
    remaining = np.ones(n, dtype=bool)
    tried = np.zeros(n, dtype=bool)
    asu_id = np.full(n, -1, dtype=int)
    num, den = as_fraction_tau(tau)

    batch_size = max(1, int(parallel_asus))
    k = 0
    while k < max_asus:
        rem_idx = np.where(remaining)[0]
        if rem_idx.size < 2:
            break

        # First filter for tracts with UR >= tau
        rem_UR = UR[rem_idx]
        high_ur_mask = rem_UR >= tau

        # If no remaining tract has UR >= tau, stop building ASUs
        if not high_ur_mask.any():
            if verbose:
                print(f"\nNo remaining tracts have UR >= {tau*100:.2f}%. Stopping ASU creation.", flush=True)
            break

        # Filter to only consider high UR tracts as potential seeds
        high_ur_rem_idx = rem_idx[high_ur_mask]

        # Among high UR tracts, find those with at least one remaining neighbor
        deg_rem = np.array([np.sum(remaining[np.array(nb[i], dtype=int)]) for i in high_ur_rem_idx])
        cand_seeds = high_ur_rem_idx[deg_rem > 0]

        if cand_seeds.size == 0:
            if verbose:
                print(f"\nNo high-UR tracts (UR >= {tau*100:.2f}%) have remaining neighbors. Stopping.", flush=True)
            break

        # Prioritize by UR (descending) then population (descending)
        order = np.lexsort((-df.loc[cand_seeds, "tract_pop2024"].to_numpy(), -UR[cand_seeds]))
        seed_pool = cand_seeds[order]

        # ---- Select up to batch_size disjoint feasible windows ----
        reserved = np.zeros(n, dtype=bool)
        windows: List[Dict] = []
        for s in seed_pool:
            if len(windows) >= batch_size:
                break
            s = int(s)
            if tried[s] or reserved[s] or not remaining[s]:
                continue
            allowed_idx = np.where(remaining & ~reserved)[0]
            if allowed_idx.size == 0:
                break

            r = int(r_start)
            sub = bfs_ball(nb, s, r, allowed_idx)
            while P[sub].sum() < min_pop_margin * pop_thresh and r < r_max and len(sub) < hard_cap_nodes:
                r += r_step
                sub = bfs_ball(nb, s, r, allowed_idx)
            if len(sub) > hard_cap_nodes:
                while len(sub) > hard_cap_nodes and r > 1:
                    r -= 1
                    sub = bfs_ball(nb, s, r, allowed_idx)
                if len(sub) > hard_cap_nodes:
                    sub = sub[:hard_cap_nodes]

            local_index = {g: i for i, g in enumerate(sub)}
            nb_local: List[List[int]] = [
                sorted(local_index[h] for h in nb[g] if h in local_index) for g in sub
            ]
            u_g, E_g, P_g = u[sub], E[sub], P[sub]
            deg_w = np.array([len(v) for v in nb_local])
            cand = np.where(deg_w > 0)[0]
            if cand.size == 0:
                tried[s] = True
                continue
            if (u_g / np.maximum(u_g + E_g, 1e-12)).max(initial=0.0) < tau:
                if verbose:
                    print(f"  [seed={s}] skip: window max(UR) < tau", flush=True)
                tried[s] = True
                continue
            if not can_hit_tau(u_g, E_g, P_g, nb_local, tau, pop_thresh):
                if verbose:
                    print(f"  [seed={s}] skip: quick screen fails", flush=True)
                tried[s] = True
                continue

            top = cand[np.argmax(u_g[cand] / np.maximum(u_g[cand] + E_g[cand], 1e-12))]
            tie = np.where(
                (u_g / np.maximum(u_g + E_g, 1e-12)) == (u_g[top] / max(u_g[top] + E_g[top], 1e-12))
            )[0]
            root_local = int(tie[np.argmax(P_g[tie])]) if len(tie) > 1 else int(top)

            windows.append({
                "seed": s, "sub": sub, "nb_local": nb_local,
                "u_g": u_g, "E_g": E_g, "P_g": P_g, "root_local": root_local, "r": r,
            })
            reserved[sub] = True

        if not windows:
            if verbose:
                print("No remaining high-UR seeds produce a feasible window; stopping.", flush=True)
            break

        if verbose:
            print(f"\n[Batch] solving {len(windows)} window(s) concurrently (ASUs {k+1}..{k+len(windows)}) ...", flush=True)

        # ---- Build warm-start hints (sequential; cheap relative to CP-SAT) ----
        for w in windows:
            info = _prepare_window_hint(
                w["nb_local"], w["u_g"], w["E_g"], w["P_g"], tau, pop_thresh, w["root_local"]
            )
            w.update(info)
            if verbose:
                su, sE = int(u[w["sub"]].sum()), int(E[w["sub"]].sum())
                URw = 100.0 * (0.0 if (su + sE) == 0 else su / (su + sE))
                seed_ur = 100.0 * UR[w["seed"]]
                print(f"\n[ASU seed={w['seed']}] (UR={seed_ur:.2f}%) | window: r={w['r']}, "
                      f"nodes={len(w['sub'])}, pop={int(w['P_g'].sum())}, UR={URw:.2f}%", flush=True)
                print(f"  root_local={w['root_local']} (UR={100*(w['u_g'][w['root_local']]/max(w['u_g'][w['root_local']]+w['E_g'][w['root_local']],1e-12)):.3f}%, "
                      f"pop={int(w['P_g'][w['root_local']])})", flush=True)
                if info["n_contracted"] < len(w["nb_local"]):
                    print(f"  contracted: {len(w['nb_local'])} -> {info['n_contracted']} nodes", flush=True)
                if len(info["root_component"]) > 1:
                    print(f"  fixed root high-UR component: {len(info['root_component'])} tracts", flush=True)
                if info["hint_valid"]:
                    hu = info["hint_obj_val"]
                    hE = int(w["E_g"][np.array(info["hint_improved"], dtype=int)].sum())
                    print(f"  [HINT] {info['hint_source']} warm start: tracts={len(info['hint_improved'])}, unemp={hu}, UR={100.0*hu/max(hu+hE,1):.2f}%", flush=True)

        # ---- Solve all windows in the batch concurrently ----
        workers_each = max(1, int(workers) // len(windows))

        def _solve(w: Dict) -> Optional[CpsatResult]:
            if verbose:
                print(f"  [seed={w['seed']}] >>> starting CP-SAT solve (nodes={len(w['nb_local'])}, "
                      f"workers={workers_each}, time_limit={time_limit}s)", flush=True)
            result = solve_one_asu_cpsat(
                nb_local=w["nb_local"], u_g=w["u_g"], E_g=w["E_g"], P_g=w["P_g"],
                tau=tau, pop_thresh=pop_thresh, root_local=w["root_local"],
                time_limit=time_limit, workers=workers_each, rel_gap=rel_gap, log=verbose,
                hint=w["hint_improved"], hint_obj=w["hint_obj_val"],
                forced_selected=w["root_component"],
                # cluster_groups intentionally NOT passed here: tying high-UR
                # cluster members via equality is provably correct (validated
                # against brute force) but empirically hurts this time-limited
                # heuristic search -- see SKILL.md "Known Issues / Gotchas".
            )
            if verbose:
                status = result.status if result is not None else "NO SOLUTION"
                print(f"  [seed={w['seed']}] <<< solve finished: status={status}", flush=True)
            return result

        if len(windows) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(windows)) as pool:
                sols = list(pool.map(_solve, windows))
        else:
            sols = [_solve(windows[0])]

        # ---- Refine each window's result within its own reserved territory ----
        batch_mask = np.zeros(n, dtype=bool)
        for w in windows:
            batch_mask[w["sub"]] = True

        candidates: List[List[int]] = []
        for w, sol in zip(windows, sols):
            if sol is None:
                if w["hint_valid"]:
                    S_local = w["hint_improved"]
                    if verbose:
                        su = int(w["u_g"][np.array(S_local, dtype=int)].sum())
                        print(f"  [seed={w['seed']}] [GREEDY FALLBACK] tracts={len(S_local)}, unemp={su}", flush=True)
                else:
                    tried[w["seed"]] = True
                    continue
            else:
                S_local = sol.sel_idx_local

            if export_dir is not None:
                _export_window_comparison(w, list(S_local), df, export_dir, k + 1)

            S_global = np.array(w["sub"], dtype=int)[np.array(S_local, dtype=int)].tolist()
            own_mask = np.zeros(n, dtype=bool)
            own_mask[w["sub"]] = True
            allowed_idx = np.where(remaining & (~batch_mask | own_mask))[0]
            S_refined = improve_by_trades(S_global, u, E, P, nb, tau, pop_thresh, allowed_idx, max_iter=200)
            if not component_ok(S_refined, u, E, P, tau, pop_thresh, nb):
                S_refined = S_global
            candidates.append(S_refined)

        if not candidates:
            continue

        # ---- Merge candidates that touch (share a queen-contiguity edge) ----
        parent = list(range(len(candidates)))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        if merge_adjacent and len(candidates) > 1:
            owner: Dict[int, int] = {}
            for gi, S in enumerate(candidates):
                for t in S:
                    owner[t] = gi
            for gi, S in enumerate(candidates):
                for t in S:
                    for w2 in nb[t]:
                        gj = owner.get(w2)
                        if gj is not None and gj != gi:
                            union(gi, gj)

        groups: Dict[int, List[int]] = {}
        for gi, S in enumerate(candidates):
            groups.setdefault(find(gi), []).append(gi)

        # Try the merged union first; fall back to each original candidate on the
        # (mathematically unexpected) chance the merged set fails a sanity check.
        final_units: List[List[int]] = []
        for members in groups.values():
            merged = sorted({t for gi in members for t in candidates[gi]})
            su, sE = int(u[merged].sum()), int(E[merged].sum())
            if len(members) > 1 and den * su - num * sE >= 0 and component_ok(merged, u, E, P, tau, pop_thresh, nb):
                if verbose:
                    print(f"  [MERGE] {len(members)} touching windows combined into one ASU ({len(merged)} tracts)", flush=True)
                final_units.append(merged)
            else:
                final_units.extend(candidates[gi] for gi in members)

        # ---- Commit each final unit, largest first, excluding not-yet-processed siblings ----
        final_units.sort(key=lambda S: -int(u[S].sum()))
        pending_mask = np.zeros(n, dtype=bool)
        for S in final_units:
            pending_mask[S] = True

        for S in final_units:
            pending_mask[S] = False
            allowed_idx = np.where(remaining & ~pending_mask)[0]
            S_final = improve_by_trades(S, u, E, P, nb, tau, pop_thresh, allowed_idx, max_iter=200)
            if not component_ok(S_final, u, E, P, tau, pop_thresh, nb):
                S_final = S

            k += 1
            asu_id[S_final] = k
            remaining[S_final] = False
            tried[S_final] = False

            if verbose:
                su, sE, sP = int(u[S_final].sum()), int(E[S_final].sum()), int(P[S_final].sum())
                URv = 100.0 * (0.0 if (su + sE) == 0 else su / (su + sE))
                print(f"  [OK] ASU {k} committed: tracts={len(S_final)}, pop={sP}, UR={URv:.3f}%, unemp={su}", flush=True)

            if k >= max_asus:
                break

    # Final summary if stopped due to no high-UR tracts
    if verbose and k < max_asus:
        rem_idx_final = np.where(remaining)[0]
        if rem_idx_final.size > 0:
            max_ur_remaining = UR[rem_idx_final].max() * 100
            print(f"\nStopped after {k} ASUs. Max UR among {rem_idx_final.size} remaining tracts: {max_ur_remaining:.3f}%", flush=True)

    return {"asu_id": asu_id, "n_asu": int(k)}


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="ASU builder with OR-Tools CP-SAT (queen contiguity supported)")
    ap.add_argument("--input", required=True, help="Excel/CSV with geoid, tract_ASU_unemp, tract_ASU_emp, tract_pop2024")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (if Excel)")
    ap.add_argument("--neighbors", default=None, help="Adjacency JSON (list of int lists; 0- or 1-based)")
    ap.add_argument("--geometry", default=None, help="GeoPackage / Shapefile with tract polygons (for queen contiguity)")
    ap.add_argument("--geom-col", default="geometry", help="Geometry column name")
    ap.add_argument("--geoid-col", default="geoid", help="Join key in geometry file (to match input geoid)")
    ap.add_argument("--tau", type=float, default=0.0645)
    ap.add_argument("--pop-thresh", type=int, default=10000)
    ap.add_argument("--max-asus", type=int, default=30)
    ap.add_argument("--r-start", type=int, default=50)
    ap.add_argument("--r-step", type=int, default=1)
    ap.add_argument("--r-max", type=int, default=50)
    ap.add_argument("--hard-cap-nodes", type=int, default=10000)
    ap.add_argument("--min-pop-margin", type=float, default=1.0)
    ap.add_argument("--time-limit", type=int, default=1200, help="CP-SAT time limit per window (seconds)")
    ap.add_argument("--workers", type=int, default=8, help="CP-SAT parallel workers")
    ap.add_argument("--rel-gap", type=float, default=None, help="Optional relative gap (e.g., 0.01 for 1%)")
    ap.add_argument("--parallel-asus", type=int, default=4, help="Number of ASU windows to solve concurrently")
    ap.add_argument("--no-merge-adjacent", action="store_true", help="Disable merging of touching ASUs built in the same batch")
    ap.add_argument("--output", default=None, help="Output CSV path (default: <stem>_with_asu.csv)")
    ap.add_argument("--verbose", action="store_true", help="Verbose CP-SAT logs")
    args = ap.parse_args()

    # Load input table
    inp = args.input
    if inp.lower().endswith((".xlsx", ".xls")):
        if args.sheet is None:
            # pick the first visible sheet
            tmp = pd.read_excel(inp, sheet_name=None)
            first_key = next(iter(tmp.keys()))
            df = tmp[first_key]
        else:
            df = pd.read_excel(inp, sheet_name=args.sheet)
    else:
        df = pd.read_csv(inp)

    # Normalize geoid (strip 14000US prefix if present)
    if "geoid" in df.columns:
        df["geoid"] = df["geoid"].astype(str).str.replace(r"^14000US", "", regex=True)

    required = ["tract_ASU_unemp", "tract_ASU_emp", "tract_pop2024"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    # Build adjacency
    if args.neighbors:
        with open(args.neighbors, "r") as f:
            nb_raw = json.load(f)
        if not isinstance(nb_raw, list):
            raise ValueError("neighbors JSON must be a list of lists")
        # Convert each row to 0-based ints; handle 1-based input from R
        n = len(nb_raw)
        nb: List[List[int]] = []
        for row in nb_raw:
            row = [int(v) for v in (row or [])]
            is_one_based = len(row) > 0 and max(row) >= n
            if is_one_based:
                row = [v - 1 for v in row]
            nb.append(sorted([v for v in row if 0 <= v < n]))
    elif args.geometry:
        if gpd is None or Queen is None:
            raise RuntimeError("geopandas/libpysal not installed. Use --neighbors JSON instead, or install geo deps.")
        gdf = gpd.read_file(args.geometry)
        if args.geoid_col not in gdf.columns:
            raise ValueError(f"Column '{args.geoid_col}' not found in geometry file.")
        # Join geometry to df by geoid
        gdf2 = gdf[[args.geoid_col, args.geom_col]].rename(columns={args.geoid_col: "geoid"})
        merged = df.merge(gdf2, on="geoid", how="left")
        if merged[args.geom_col].isna().any():
            missing = merged["geoid"][merged[args.geom_col].isna()].unique()[:5]
            raise RuntimeError(f"Missing geometry for some geoids (e.g., {missing}).")
        gdf_merged = gpd.GeoDataFrame(merged, geometry=args.geom_col, crs=gdf.crs).reset_index(drop=True)
        nb = queen_neighbors_from_geometries(gdf_merged, geom_col=args.geom_col)
        # Drop geometry for output size
        df = pd.DataFrame(gdf_merged.drop(columns=[args.geom_col]))
    else:
        raise RuntimeError("Provide either --neighbors JSON or --geometry to compute contiguity.")

    # Build ASUs
    out = build_many_asus_cpsat(
        df=df, nb=nb, tau=args.tau, pop_thresh=args.pop_thresh,
        max_asus=args.max_asus, r_start=args.r_start, r_step=args.r_step, r_max=args.r_max,
        hard_cap_nodes=args.hard_cap_nodes, min_pop_margin=args.min_pop_margin,
        time_limit=args.time_limit, workers=args.workers, rel_gap=args.rel_gap,
        verbose=args.verbose, parallel_asus=args.parallel_asus,
        merge_adjacent=not args.no_merge_adjacent,
    )

    df_out = df.copy()
    df_out["asu_id"] = out["asu_id"]

    out_path = args.output or f"{os.path.splitext(os.path.basename(inp))[0]}_with_asu.csv"
    df_out.to_csv(out_path, index=False)
    print(f"\nDone. Built {out['n_asu']} ASU(s) → {out_path}")


if __name__ == "__main__":
    main()
