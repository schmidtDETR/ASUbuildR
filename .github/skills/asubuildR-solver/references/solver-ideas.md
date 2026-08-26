# Solver Improvement Ideas — Detailed Notes

## Validated Changes (2026-07-30)

- Root high-UR component members are fixed selected using the exact integer UR
    inequality. Exhaustive randomized checks confirmed the optimal objective is preserved.
- A bounded separator-cut phase precedes the exact flow model. On Colorado 2025 at
    300 seconds: objective 74,805, bound 75,690, versus the prior 74,088 / 75,818.
- `improve_by_trades` now maintains incremental unemployment, employment, and population
    totals. A 60-second end-to-end benchmark fell from 230.1 seconds to 63.3 seconds.
- Hybrid model exactness was checked against exhaustive enumeration on 60 randomized
    connected graphs; the incremental improver matched its predecessor on 300 graphs.

---

## 1. Fix x-variables for Root Super-Node Members

**What**: After `contract_high_ur_nodes`, `expand_r[root_r]` contains all original local
indices of high-UR tracts connected to the root. These tracts MUST all be selected in any
feasible ASU containing root (mediant inequality: mixing two UR≥tau groups preserves UR≥tau,
so they're inseparable from root's cluster). Fix them in the model:

```python
root_supernode_members = expand_r[root_r]  # list of original local indices
for v in root_supernode_members:
    model.Add(x[v] == 1)
```

**Why it helps**: Fixes variables → presolve propagation eliminates them → smaller LP,
tighter flow bounds, fewer branch-and-bound nodes.

**Risk**: None. If the super-node is a singleton (root is alone), the list is `[root_local]`
which is already forced by `model.Add(x[root_local] == 1)`.

**Implementation note**: These fixed variables still need to appear in the hint to avoid
"complete hint is infeasible" warning from CP-SAT. They already do (all members of the root
super-node are in `hint_expanded`).

---

## 2. Tighter Per-Edge Flow Bounds

**What**: Currently every flow variable uses `M = N-1`. Instead compute subtree sizes from
a BFS spanning tree rooted at `root_local` and use `subtree_size[j]` as the upper bound
for directed edge `(i→j)`.

```python
# Compute BFS spanning tree from root_local (using full nb_local)
subtree_size = {v: 1 for v in range(N)}
parent = {}
bfs_order = []
queue = [root_local]; visited = {root_local}
children = {v: [] for v in range(N)}
while queue:
    v = queue.pop(0); bfs_order.append(v)
    for w in nb_local[v]:
        if w not in visited:
            visited.add(w); parent[w] = v; children[v].append(w); queue.append(w)
for v in reversed(bfs_order):
    if v in parent: subtree_size[parent[v]] += subtree_size[v]

# Build flow variables with tight bounds
f = []
for (i, j) in edges:
    # Flow on (i→j) is bounded by: how many nodes are "behind j from i"
    # In the BFS tree, if j is a child of i: bound = subtree_size[j]
    # If i is a child of j (back-edge): bound = N - subtree_size[i]
    # Otherwise (cross-edge): bound = N-1 (conservative)
    if parent.get(j) == i:
        ub = subtree_size[j]
    elif parent.get(i) == j:
        ub = N - subtree_size[i]
    else:
        ub = N - 1
    f.append(model.NewIntVar(0, max(1, ub), f"f_{i}_{j}"))
```

**Why it helps**: Tighter M → tighter LP relaxation → faster bound convergence.
Tree edges (majority of edges in sparse planar graphs) get bounds of O(subtree_size)
instead of O(N). For a balanced tree this halves the average bound.

---

## 3. Solve on Contracted Graph with Equality Links

**What**: Build the model on contracted variables `xr[ri]` for `ri in range(N_r)`,
then link to original with `x[orig] == xr[node_map[orig]]` for all orig.
Effectively, the solver only branches on N_r decision variables instead of N.

```python
xr = [model.NewBoolVar(f"xr_{ri}") for ri in range(N_r)]
x  = [model.NewBoolVar(f"x_{v}")   for v in range(N)]
for v in range(N):
    model.Add(x[v] == xr[int(node_map[v])])
# Use xr in all constraints; x is auxiliary for hint passing
```

**Why it helps**: For states where many tracts are UR≥tau and cluster together,
N_r << N. Branching on xr only → dramatically fewer nodes in B&B tree.

**Pitfall**: Flow variables must be defined over the CONTRACTED adjacency `nb_r`,
not the original. The flows on contracted edges represent bundled connectivity.
Implement carefully — verify on CO test first.

---

## 4. Adaptive Per-ASU Time Limit

**What**: Scale time limit by relative window size:

```python
frac = min(1.0, len(sub) / max(1, n_total))
effective_limit = max(30, int(time_limit * (frac ** 0.7)))
```

**Why it helps**: For a 200-node window in a 1447-node state, there's no reason to
spend 300s. Freeing time for genuinely hard large windows matters in multi-ASU states.

---

## 5. Two-Phase Solve

**What**: Run CP-SAT for 20% of time_limit first, record the incumbent objective,
then restart with `hint_obj = incumbent + delta` and remaining 80% of time.

Actually simpler: CP-SAT already does this internally. The real gain is if the greedy
hint is known to be far from optimal (e.g., hint achieves 50% of LP bound), in which
case setting `rel_gap = 0.02` (2%) terminates early when within 2% of LP bound.

```python
solver.parameters.relative_gap_limit = 0.02  # accept 2% suboptimality
```

**Tradeoff**: Relaxing rel_gap slightly reduces solution quality but can save 50-80%
of solver time on hard instances. Currently rel_gap=None (find exact optimum within
time_limit). Consider making it a configurable parameter.

---

## 6. Presolve: Fix Zero-Contribution Tracts

**What**: Before building the model, identify tracts that can NEVER contribute to
feasibility or improve the objective, and fix `x[v] = 0`:

```python
# Tracts with zero unemployment never improve objective
zero_u = np.where(u_g == 0)[0]
for v in zero_u:
    model.Add(x[v] == 0)
```

More aggressively: any tract v where `u_g[v] == 0 AND P_g[v] == 0` is pure dead weight.
Any tract v that is isolated in the window (no neighbors) cannot be reached by flow → fix 0.

**Why it helps**: Fewer active variables → faster presolve propagation.

---

## 7. Parallel Seed Evaluation

**What**: When the first seed fails quick screen, evaluate the next K seeds'
`can_hit_tau` in parallel (Python `concurrent.futures.ThreadPoolExecutor`) before
entering the serial solve loop.

```python
from concurrent.futures import ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=4) as pool:
    futures = {pool.submit(can_hit_tau, u_g, E_g, P_g, nb_local, tau, pop_thresh): seed
               for seed, (u_g, E_g, P_g, nb_local) in candidate_windows.items()}
    valid_seeds = [seed for f, seed in futures.items() if f.result()]
```

**Note**: `can_hit_tau` is pure Python/NumPy with no shared state — safe to parallelize.

---

## 8. Post-ASU Look-Ahead

**What**: After committing ASU k, check all candidate seeds for ASU k+1 BEFORE
running the solver. If NO seed passes `can_hit_tau`, skip the solve loop entirely.

Currently the code does check `can_hit_tau` per seed but only when that seed is
popped from the queue. A bulk pre-check over all remaining high-UR seeds avoids
entering the main loop N_seeds times only to fail each time.

```python
# After committing ASU k
rem_high = np.where(remaining & (UR >= tau))[0]
viable = any(
    can_hit_tau(u[bfs_ball(nb, s, r_start, rem_high)], ...)
    for s in rem_high[:20]  # check top 20 seeds
)
if not viable:
    break  # skip solve loop
```

---

## CP-SAT Model Reference

```
Variables: x[i] ∈ {0,1},  f[e] ∈ [0, M]  (e = directed edge (i,j))
M = N-1  (current; improve with per-edge subtree bounds — idea #2)

Constraints:
  x[root] = 1
  Σ P[i]*x[i] ≥ pop_thresh
  den*Σ u[i]*x[i] - num*Σ E[i]*x[i] ≥ 0      (UR ≥ tau, integer exact)
  f[e] ≤ M*x[i],  f[e] ≤ M*x[j]              (flow only on selected edges)
  For root:   Σ_out f - Σ_in f = Σ x[i] - 1  (source)
  For i≠root: Σ_in f - Σ_out f = x[i]         (sink if selected)

Objective: max Σ u[i]*x[i]

Hints: AddHint(x[i], hint_bit) + AddHint(f[e], spanning_tree_flow)
Lower bound: model.Add(Σ u[i]*x[i] ≥ hint_obj)  [rejects worse solutions]

Solver params:
  linearization_level = 2
  cp_model_presolve = True
  log_search_progress = verbose
  num_search_workers = workers (default 18)
  max_time_in_seconds = time_limit (default 300)
```

---

## Data Schema

Input DataFrame columns:
- `geoid` — Census tract GEOID (string)
- `tract_ASU_unemp` — unemployed persons (int)
- `tract_ASU_emp` — employed persons (int)
- `tract_pop2024` — total population (int, used for pop_thresh check)

`nb` — list of lists (0-based), length = number of tracts. `nb[i]` = list of
0-based indices of queen-contiguous neighbors of tract i. Precomputed in R via
`sfdep::st_contiguity()` and written as JSON.

---

## Function Index

| Function | File | Purpose |
|----------|------|---------|
| `build_many_asus_cpsat` | asu_cpsat.py L560+ | Top-level loop: seed → window → solve → commit |
| `solve_one_asu_cpsat` | asu_cpsat.py L360+ | Build and solve single CP-SAT model |
| `greedy_snake_hint` | asu_cpsat.py L84+ | Warm-start hint via Snake + group merge |
| `improve_by_trades` | asu_cpsat.py L510+ | Post-solve greedy add + swap improver |
| `contract_high_ur_nodes` | asu_cpsat.py L255+ | Fuse UR≥tau clusters into super-nodes |
| `can_hit_tau` | asu_cpsat.py L240+ | Fractional knapsack feasibility screen |
| `component_ok` | asu_cpsat.py L228+ | Connectivity + UR/pop check |
| `_spanning_tree_flows` | asu_cpsat.py L195+ | Flow hints from BFS spanning tree |
| `bfs_ball` | asu_cpsat.py L60+ | BFS window around seed |
| `as_fraction_tau` | asu_cpsat.py L55+ | tau → (num, den) exact integer fraction |
