#!/usr/bin/env python3

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"


def run_cpp_alg(alg_name: str, budget: int, rows: int = None, cols: int = None, seed: int = None,
                input_csv_path: Path = None):
    if not CPP_BIN.exists():
        raise RuntimeError(f"C++ binary not found at {CPP_BIN}")

    cmd = [
        str(CPP_BIN),
        "-algorithm", alg_name,
        "-budget", str(budget),
        "-log", "0",
    ]
    if input_csv_path is not None:
        cmd.extend(["-input_csv", str(input_csv_path)])
    else:
        if rows is None or cols is None or seed is None:
            raise ValueError("rows/cols/seed are required when input_csv_path is not provided")
        cmd.extend(["-rows", str(rows), "-cols", str(cols), "-seed", str(seed)])
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, capture_output=True, text=True)

    reward = None
    cost = None
    cycle = []
    partial_row_feasible = None
    partial_row_reason = ""
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            reward = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            cost = float(line.split("=", 1)[1].strip())
        elif line.startswith("partial_row_feasible="):
            partial_row_feasible = int(line.split("=", 1)[1].strip())
        elif line.startswith("partial_row_reason="):
            partial_row_reason = line.split("=", 1)[1].strip()
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
        "partial_row_feasible": partial_row_feasible,
        "partial_row_reason": partial_row_reason,
    }


def draw_instance_with_path(ax, grid, cycle, title, reward_value, cost_value,
                            partial_row_feasible=None, partial_row_reason=""):
    internal = np.array(grid, dtype=float)
    m, n = internal.shape

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
    if partial_row_feasible is not None:
        status = "OK" if int(partial_row_feasible) == 1 else "KO"
        if partial_row_reason:
            status += f" ({partial_row_reason})"
    title_lines = [f"{title}", f"reward={reward_value:.2f}, cost={cost_value:.2f}"]
    if status:
        title_lines.append(status)
    ax.set_title("\n".join(title_lines), fontsize=10,
                 color=("#1b7f3b" if status.startswith("OK") else ("#b42318" if status else "black")))
    ax.set_xlabel("graph column (0..n+1)")
    ax.set_ylabel("row i")
    ax.set_xticks(range(n + 2))
    ax.set_yticks(range(m))
    ax.set_xlim(-0.5, n + 1.5)
    ax.set_ylim(m - 0.5, -0.5)
    ax.grid(color="white", linestyle="-", linewidth=0.8, alpha=0.5)


def parse_args():
    p = argparse.ArgumentParser(description="Plot OPRSC vs GPRSC paths on one aisle-grid instance (internal columns only).")
    p.add_argument("--budget", type=int, required=True, help="Travel budget B")
    p.add_argument("--input-csv", type=str, default="", help="CSV with m x n internal rewards")
    p.add_argument("--rows", type=int, default=6, help="Rows m (used if --input-csv is omitted)")
    p.add_argument("--cols", type=int, default=8, help="Internal columns n (used if --input-csv is omitted)")
    p.add_argument("--seed", type=int, default=0, help="Random seed (used if --input-csv is omitted)")
    p.add_argument("--out", type=str, default="test/output/plot_oprsc_gprsc.pdf", help="Output figure path (e.g., .pdf, .png)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.input_csv:
        with open(args.input_csv, "r", newline="") as f:
            grid = [[float(x) for x in row] for row in csv.reader(f) if row]
    else:
        rng = np.random.default_rng(args.seed)
        grid = rng.integers(0, 11, size=(args.rows, args.cols)).astype(float).tolist()

    if not grid or not grid[0]:
        raise ValueError("Empty grid")

    m = len(grid)
    n = len(grid[0])
    if any(len(row) != n for row in grid):
        raise ValueError("Grid must be rectangular")

    if args.input_csv:
        input_csv = Path(args.input_csv)
        oprsc_cpp = run_cpp_alg("oprsc", args.budget, input_csv_path=input_csv)
        gprsc_cpp = run_cpp_alg("gprsc", args.budget, input_csv_path=input_csv)
    else:
        oprsc_cpp = run_cpp_alg("oprsc", args.budget, rows=args.rows, cols=args.cols, seed=args.seed)
        gprsc_cpp = run_cpp_alg("gprsc", args.budget, rows=args.rows, cols=args.cols, seed=args.seed)

    if oprsc_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"OPRSC infeasible: cost={oprsc_cpp['cost']} > budget={args.budget}")
    if gprsc_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"GPRSC infeasible: cost={gprsc_cpp['cost']} > budget={args.budget}")
    if oprsc_cpp.get("partial_row_feasible") != 1:
        raise RuntimeError(f"OPRSC validator KO: {oprsc_cpp.get('partial_row_reason', '')}")
    if gprsc_cpp.get("partial_row_feasible") != 1:
        raise RuntimeError(f"GPRSC validator KO: {gprsc_cpp.get('partial_row_reason', '')}")
    if gprsc_cpp["reward"] > oprsc_cpp["reward"] + 1e-9:
        raise RuntimeError(
            "Ordering violated: GPRSC reward exceeds OPRSC reward "
            f"(gprsc={gprsc_cpp['reward']}, oprsc={oprsc_cpp['reward']}, budget={args.budget})"
        )

    fig, axes = plt.subplots(1, 2, figsize=(max(10, n * 1.6), max(4.8, m * 0.8)), constrained_layout=True)
    fig.suptitle(
        f"Single-instance OPRSC vs GPRSC (m={m}, n={n} internal cols, graph cols={n + 2}, B={args.budget})\n"
        "Note: corridor columns are implicit and have zero profit",
        fontsize=12,
    )

    draw_instance_with_path(
        axes[0], grid, oprsc_cpp["cycle"], "OPRSC", float(oprsc_cpp["reward"]), float(oprsc_cpp["cost"]),
        oprsc_cpp.get("partial_row_feasible"), oprsc_cpp.get("partial_row_reason", "")
    )
    draw_instance_with_path(
        axes[1], grid, gprsc_cpp["cycle"], "GPRSC", float(gprsc_cpp["reward"]), float(gprsc_cpp["cost"]),
        gprsc_cpp.get("partial_row_feasible"), gprsc_cpp.get("partial_row_reason", "")
    )

    out_path = Path(args.out)
    if not out_path.parent.exists():
        raise RuntimeError(f"Output directory does not exist: {out_path.parent}")
    if not out_path.parent.is_dir():
        raise RuntimeError(f"Output path parent is not a directory: {out_path.parent}")
    fig.savefig(out_path, dpi=180)
    print(f"Saved figure to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
