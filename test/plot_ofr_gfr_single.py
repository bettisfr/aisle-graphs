#!/usr/bin/env python3

import argparse
import csv
import math
import random
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
OLD_INPUT_DIR = REPO_ROOT / "_old" / "input"
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"


def write_grid_csv(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        csv.writer(f).writerows(grid)


def run_cpp_alg(alg_name: str, input_csv_path: Path, budget: int):
    if not CPP_BIN.exists():
        raise RuntimeError(f"C++ binary not found at {CPP_BIN}")

    cmd = [
        str(CPP_BIN),
        "-algorithm", alg_name,
        "-budget", str(budget),
        "-input_csv", str(input_csv_path),
        "-log", "0",
    ]
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, capture_output=True, text=True)

    reward = None
    cost = None
    cycle = []
    full_row_feasible = None
    full_row_reason = ""
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            reward = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            cost = float(line.split("=", 1)[1].strip())
        elif line.startswith("full_row_feasible="):
            full_row_feasible = int(line.split("=", 1)[1].strip())
        elif line.startswith("full_row_reason="):
            full_row_reason = line.split("=", 1)[1].strip()
        elif line.startswith("cycle_list="):
            raw = line.split("=", 1)[1].strip()
            if raw:
                for token in raw.split(";"):
                    if not token:
                        continue
                    r_str, c_str = token.split(":")
                    cycle.append((int(r_str), int(c_str)))
    if reward is None or cost is None:
        raise RuntimeError(f"Cannot parse {alg_name} output.\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

    return {
        "reward": reward,
        "cost": cost,
        "cycle": cycle,
        "full_row_feasible": full_row_feasible,
        "full_row_reason": full_row_reason,
    }


def build_random_grid(rows: int, cols: int, seed: int, min_reward: int = 0, max_reward: int = 10):
    rng = random.Random(seed)
    return [[rng.randint(min_reward, max_reward) for _ in range(cols)] for _ in range(rows)]


def greedy_full_row_reference_with_path(grid, budget: int):
    r = np.array(grid, dtype=float)
    m, n = r.shape

    out = {"reward": 0.0, "cost": 0.0, "traversed_rows": [], "cycle": []}
    if budget <= 0 or m == 0 or n == 0:
        out["cycle"] = [(0, 0)]
        return out

    compiled = [(float(np.sum(r[i, :])), i + 1) for i in range(m)]  # 1-based rows
    compiled.sort(key=lambda x: (-x[0], x[1]))

    beginning_row = 1
    ending_row = 1
    compiled_row_cost = n
    current_row = beginning_row
    current_side = 1
    total_cost = 0
    total_reward = 0.0
    selected_rows_1based = []

    for row_reward, row_id in compiled:
        if current_side == 1:
            loop_cost = abs(current_row - row_id) + 2 * compiled_row_cost + abs(row_id - ending_row)
        else:
            loop_cost = abs(current_row - row_id) + compiled_row_cost + abs(row_id - ending_row)

        if loop_cost <= (budget - total_cost):
            total_cost += abs(current_row - row_id) + compiled_row_cost
            total_reward += row_reward
            selected_rows_1based.append(row_id)
            current_row = row_id
            current_side = 2 if current_side == 1 else 1

    if current_side == 2:
        total_cost += abs(current_row - ending_row) + compiled_row_cost
    else:
        total_cost += abs(current_row - ending_row)

    out["reward"] = float(total_reward)
    out["cost"] = float(total_cost)
    out["traversed_rows"] = [rid - 1 for rid in selected_rows_1based]

    # Compressed cycle on internal columns (same convention as C++ port)
    left_col = 0
    right_col = n - 1
    side_col = left_col
    cycle = [(0, side_col)]
    for rid in selected_rows_1based:
        rr = rid - 1
        if cycle[-1] != (rr, side_col):
            cycle.append((rr, side_col))
        side_col = right_col if side_col == left_col else left_col
        if cycle[-1] != (rr, side_col):
            cycle.append((rr, side_col))
    if cycle[-1] != (0, side_col):
        cycle.append((0, side_col))
    if side_col != left_col:
        cycle.append((0, left_col))

    out["cycle"] = cycle
    return out


def draw_instance_with_path(ax, grid, cycle, title, reward_value, cost_value, full_row_feasible=None, full_row_reason=""):
    internal = np.array(grid, dtype=float)
    m, n = internal.shape

    # Build full aisle-graph reward view: corridor columns c_0 and c_{n+1} have reward 0.
    full = np.zeros((m, n + 2), dtype=float)
    full[:, 1:n + 1] = internal

    ax.imshow(
        full,
        cmap="Pastel1",
        interpolation="nearest",
        vmin=np.min(full),
        vmax=np.max(full) if np.max(full) > 0 else 1,
        alpha=0.55,
    )

    # Draw node values on all columns, including explicit zero-profit corridors.
    for i in range(m):
        for j in range(n + 2):
            val = full[i, j]
            ax.text(j, i, f"{int(val)}", ha="center", va="center", fontsize=10, color="#0b4f9c", fontweight="bold")

    if cycle:
        xs = [c for _, c in cycle]
        ys = [r for r, _ in cycle]
        ax.plot(xs, ys, color="#f4a6a6", linewidth=2.6, marker="o", markersize=3.8, alpha=0.95)
        ax.scatter([xs[0]], [ys[0]], color="#1f77b4", s=45, zorder=5, label="start")

    status = ""
    if full_row_feasible is not None:
        status = "OK" if int(full_row_feasible) == 1 else "KO"
        if full_row_reason:
            status += f" ({full_row_reason})"
    title_lines = [f"{title}", f"reward={reward_value:.2f}, cost={cost_value:.2f}"]
    if status:
        title_lines.append(status)
    ax.set_title("\n".join(title_lines), fontsize=10, color=("#1b7f3b" if status.startswith("OK") else ("#b42318" if status else "black")))
    ax.set_xlabel("graph column (0..n+1)")
    ax.set_ylabel("row i")
    ax.set_xticks(range(n + 2))
    ax.set_yticks(range(m))
    ax.set_xlim(-0.5, n + 1.5)
    ax.set_ylim(m - 0.5, -0.5)
    ax.grid(color="white", linestyle="-", linewidth=0.8, alpha=0.5)


def parse_args():
    p = argparse.ArgumentParser(description="Plot OFR vs GFR paths on one aisle-grid instance (internal columns only).")
    p.add_argument("--budget", type=int, required=True, help="Travel budget B")
    p.add_argument("--input-csv", type=str, default="", help="CSV with m x n internal rewards")
    p.add_argument("--rows", type=int, default=6, help="Rows m (used if --input-csv is omitted)")
    p.add_argument("--cols", type=int, default=8, help="Internal columns n (used if --input-csv is omitted)")
    p.add_argument("--seed", type=int, default=0, help="Random seed (used if --input-csv is omitted)")
    p.add_argument("--out", type=str, default="test/output/plot_ofr_gfr_single.pdf", help="Output figure path (e.g., .pdf, .png)")
    p.add_argument("--legacy-input-id", type=int, default=990010, help="Temporary legacy input id for OFR reference")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.input_csv:
        with open(args.input_csv, "r", newline="") as f:
            grid = [[float(x) for x in row] for row in csv.reader(f) if row]
    else:
        grid = build_random_grid(args.rows, args.cols, args.seed)

    if not grid or not grid[0]:
        raise ValueError("Empty grid")

    m = len(grid)
    n = len(grid[0])
    if any(len(row) != n for row in grid):
        raise ValueError("Grid must be rectangular")

    # Save shared instance to CSV consumed by the C++ binary.
    legacy_input_file = OLD_INPUT_DIR / f"input-{args.legacy_input_id}.csv"
    write_grid_csv(legacy_input_file, grid)

    ofr_cpp = run_cpp_alg("ofr", legacy_input_file, args.budget)
    gfr_cpp = run_cpp_alg("gfr", legacy_input_file, args.budget)

    if ofr_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"OFR infeasible: cost={ofr_cpp['cost']} > budget={args.budget}")
    if gfr_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"GFR infeasible: cost={gfr_cpp['cost']} > budget={args.budget}")
    if gfr_cpp["reward"] > ofr_cpp["reward"] + 1e-9:
        raise RuntimeError(
            "Ordering violated: GFR reward exceeds OFR reward "
            f"(gfr={gfr_cpp['reward']}, ofr={ofr_cpp['reward']}, budget={args.budget})"
        )

    fig, axes = plt.subplots(1, 2, figsize=(max(10, n * 1.6), max(4.8, m * 0.8)), constrained_layout=True)
    fig.suptitle(
        f"Single-instance OFR vs GFR (m={m}, n={n} internal cols, graph cols={n + 2}, B={args.budget})\n"
        "Note: corridor columns are implicit and have zero profit",
        fontsize=12,
    )

    draw_instance_with_path(
        axes[0], grid, ofr_cpp["cycle"], "OFR", float(ofr_cpp["reward"]), float(ofr_cpp["cost"]),
        ofr_cpp.get("full_row_feasible"), ofr_cpp.get("full_row_reason", "")
    )
    draw_instance_with_path(
        axes[1], grid, gfr_cpp["cycle"], "GFR", float(gfr_cpp["reward"]), float(gfr_cpp["cost"]),
        gfr_cpp.get("full_row_feasible"), gfr_cpp.get("full_row_reason", "")
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    print(f"Saved figure to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
