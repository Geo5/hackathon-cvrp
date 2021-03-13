import datetime as dt
import json
import sys

import numpy as np
import scipy.sparse as sparse
from lpsolve55 import *
from scipy.sparse.csgraph import shortest_path


def load_instance(json_file):
    inst = json.load(open(json_file))
    n = inst["graph"]["nodes"]
    edges = np.array(inst["graph"]["edges"])
    graph = np.zeros((n, n), dtype=float)
    for e in edges:
        graph[e[0] - 1, e[1] - 1] = float(e[2])
        graph[e[1] - 1, e[0] - 1] = float(e[2])
    distances = sparse.csr_matrix(graph)
    distances = shortest_path(distances, directed=False)
    pickpool = inst["pickpool"]
    # Duplicate every node which has more than 1 pickpool entry
    pick_count = {}
    mod_n = n
    for i, p in enumerate(pickpool):
        if p[0] - 1 not in pick_count:
            pick_count[p[0] - 1] = []
        else:
            mod_n += 1
        pick_count[p[0] - 1].append((p[1], i))
    mod_dists = np.zeros((mod_n, mod_n), dtype=float)
    mod_dists[0:n, 0:n] = distances
    added_rows = 0
    added_row_idxs = []
    mod_pickpool = np.zeros((mod_n,), dtype=float)
    pickpool_orig_idx = np.zeros((len(pickpool),), dtype=int)
    for i, (node, weights_idx) in enumerate(pick_count.items()):
        mod_pickpool[node] = float(weights_idx[0][0])
        pickpool_orig_idx[node - 1] = weights_idx[0][1]
        for j, w_i in enumerate(weights_idx[1:]):
            mod_dists[n : n + added_rows, n + added_rows] = distances[
                added_row_idxs, node
            ]
            mod_dists[n + added_rows, n : n + added_rows] = distances[
                added_row_idxs, node
            ]
            mod_dists[n + added_rows, :n] = distances[node, :]
            mod_dists[:n, n + added_rows] = distances[:, node]
            added_row_idxs.append(node)
            mod_pickpool[n + added_rows] = float(w_i[0])
            pickpool_orig_idx[len(pick_count.items()) + j] = w_i[1]
            added_rows += 1
    # Remove every node, that has 0 pickpool entries
    delete = mod_pickpool < 1
    delete[0] = False
    mod_dists = np.delete(mod_dists, delete, axis=0)
    mod_dists = np.delete(mod_dists, delete, axis=1)
    demands = np.delete(mod_pickpool, delete)[1:]
    orig_indices = np.delete(np.array(list(range(n)) + added_row_idxs), delete)
    return mod_dists, orig_indices, pickpool_orig_idx, demands, float(inst["cap"])


def compute_savings(distances):
    n = distances.shape[0]
    savings = np.full_like(distances, -1)
    for i in range(1, n):
        for j in range(1, n):
            if i == j:
                continue
            s = distances[0, i] + distances[0, j] - distances[i, j]
            savings[i, j] = s
            savings[j, i] = s
    return savings


def main(args):
    start_time = dt.datetime.now()
    distances, orig_indices, pickpool_orig_idx, demands, cap = load_instance(args[0])
    p = np.ceil(sum(demands) / cap)
    n = distances.shape[0]
    savings = compute_savings(distances)

    col_names = [
        "x" + str(i) + "j" + str(j) for i in range(1, n + 1) for j in range(1, n + 1)
    ]
    col_names.extend("y" + str(i) for i in range(2, n + 1))
    col_names = np.array(col_names, dtype=str)
    n_cols = len(col_names)

    lp = lpsolve("make_lp", 0, n_cols)
    for i in range(n_cols):
        lpsolve("set_col_name", lp, i + 1, col_names[i])

    obj_values = savings.flatten(order="C").clip(min=0).tolist() + [0.0] * (n - 1)
    lpsolve("set_obj_fn", lp, obj_values)
    lpsolve("set_maxim", lp)

    # Add constraints
    lpsolve("set_add_rowmode", lp, True)
    fst_non_zero = list(range(1, n))
    fst_coeff = np.zeros((n_cols,), dtype=float)
    fst_coeff[fst_non_zero] = 1
    lpsolve("add_constraintex", lp, fst_coeff, "EQ", p)
    lpsolve("set_row_name", lp, 1, "FIRST")

    for j in range(1, n):
        snd_non_zero = [i for i in range(j, n * n, n) if i // n != j]
        snd_count = len(snd_non_zero)
        snd_coeff = np.zeros((n_cols,), dtype=float)
        snd_coeff[snd_non_zero] = 1
        lpsolve("add_constraintex", lp, snd_coeff, "EQ", 1.0)
        lpsolve("set_row_name", lp, 1 + j, "SECOND" + str(j + 1))
    for i in range(1, n):
        thr_non_zero = [j for j in range(i * n + 1, (i + 1) * n) if i != j % n]
        thr_coeff = np.zeros((n_cols,), dtype=float)
        thr_coeff[thr_non_zero] = 1
        lpsolve("add_constraintex", lp, thr_coeff, "LE", 1.0)
        lpsolve("set_row_name", lp, n + i, "THIRD" + str(i + 1))
    for i in range(1, n):
        for j in range(1, n):
            if i == j:
                continue
            fth_non_zero = [i * n + j, n * n + min(i, j) - 1, n * n + max(i, j) - 1]
            fth_coeff = np.zeros((n_cols,), dtype=float)
            if i < j:
                fth_coeff[fth_non_zero] = [demands[j - 1] + cap, 1.0, -1.0]
            else:
                fth_coeff[fth_non_zero] = [demands[j - 1] + cap, -1.0, 1.0]
            lpsolve(
                "add_constraintex",
                lp,
                fth_coeff,
                "LE",
                cap,
            )
            lpsolve(
                "set_row_name",
                lp,
                n + n - 1 + (i - 1) * (n - 2) + j - int(j > i),
                "FOURTH" + str(i + 1) + str(j + 1),
            )
    for i in range(1, n):
        lpsolve("set_bounds", lp, n * n + i, demands[i - 1], cap)
    for i in range(n * n):
        lpsolve("set_binary", lp, i + 1, True)

    # Writing and re-reading of the lp, because of a bug in lpsolve(?)
    lpsolve("write_lp", lp, "test.lp")
    lpsolve("delete_lp", lp)

    lp = lpsolve("read_lp", "test.lp")
    lpsolve("set_verbose", lp, "CRITICAL")
    time_left = 100 - (dt.datetime.now() - start_time).total_seconds()
    # Set timeout, so the programm will return its best solution before 100 seconds.
    lpsolve("set_timeout", lp, max(1, time_left - 2))

    ret = lpsolve("solve", lp)
    if ret not in [0, 1]:
        print("[]")
        return

    max_savings = lpsolve("get_objective", lp)
    min_length = distances[0, 1:].sum() * 2 - max_savings

    vars = lpsolve("get_variables", lp)[0]
    edges = []
    for i_v in range(len(vars)):
        if vars[i_v] > 1e-6:
            name = lpsolve("get_col_name", lp, i_v + 1)
            if name.startswith("x"):
                i, j = name[1:].split("j")
                edges.append((int(i) - 1, int(j) - 1))
    tour_starts = [e for e in edges if e[0] == 0]
    tours = []
    for i_t in range(len(tour_starts)):
        cur_e = tour_starts[i_t]
        tour = [pickpool_orig_idx[cur_e[1] - 1]]
        cur_e = next((e for e in edges if e[0] == cur_e[1]), None)
        while cur_e is not None:
            tour.append(pickpool_orig_idx[cur_e[1] - 1])
            cur_e = next((e for e in edges if e[0] == cur_e[1]), None)
        tours.append(tour)

    print(tours)
    lpsolve("delete_lp", lp)


if __name__ == "__main__":
    main(sys.argv[1:])
