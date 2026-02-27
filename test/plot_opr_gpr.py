#!/usr/bin/env python3

import argparse
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

    out = {
        "reward": None,
        "cost": None,
        "cycle": [],
        "grid": [],
        "partial_row_feasible": None,
        "partial_row_reason": "",
    }
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            out["reward"] = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            out["cost"] = float(line.split("=", 1)[1].strip())
        elif line.startswith("partial_row_feasible="):
            out["partial_row_feasible"] = int(line.split("=", 1)[1].strip())
        elif line.startswith("partial_row_reason="):
            out["partial_row_reason"] = line.split("=", 1)[1].strip()
        elif line.startswith("cycle_list="):
            raw = line.split("=", 1)[1].strip()
            if raw:
                for token in raw.split(";"):
                    if not token:
                        continue
                    r_str, c_str = token.split(":")
                    out["cycle"].append((int(r_str), int(c_str)))
        elif line.startswith("grid_row="):
            raw = line.split("=", 1)[1].strip()
            out["grid"].append([float(x) for x in raw.split(",") if x != ""])

    if out["reward"] is None or out["cost"] is None:
        raise RuntimeError(f"Cannot parse {alg_name} output.\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return out


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
        ax.scatter([xs[0]], [ys[0]], color="#1f77b4", s=45, zorder=5)

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
    p = argparse.ArgumentParser(description="Plot OPR vs GPR paths on one aisle-grid instance.")
    p.add_argument("--budget", type=int, required=True, help="Travel budget B")
    p.add_argument("--input-csv", type=str, default="", help="CSV with m x n internal rewards")
    p.add_argument("--rows", type=int, default=10, help="Rows m (used if --input-csv is omitted)")
    p.add_argument("--cols", type=int, default=10, help="Internal columns n (used if --input-csv is omitted)")
    p.add_argument("--seed", type=int, default=0, help="Random seed (used if --input-csv is omitted)")
    p.add_argument("--out", type=str, default="output/plot_opr_gpr.pdf", help="Output figure path")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.input_csv:
        input_csv = Path(args.input_csv)
        opr_cpp = run_cpp_alg("opr", args.budget, input_csv_path=input_csv)
        gpr_cpp = run_cpp_alg("gpr", args.budget, input_csv_path=input_csv)
    else:
        opr_cpp = run_cpp_alg("opr", args.budget, rows=args.rows, cols=args.cols, seed=args.seed)
        gpr_cpp = run_cpp_alg("gpr", args.budget, rows=args.rows, cols=args.cols, seed=args.seed)

    grid = opr_cpp.get("grid", [])
    if not grid:
        raise RuntimeError("Cannot parse grid from C++ output")
    if gpr_cpp.get("grid") and gpr_cpp["grid"] != grid:
        raise RuntimeError("C++ returned different random grids across OPR/GPR calls")

    m = len(grid)
    n = len(grid[0]) if m > 0 else 0
    if m == 0 or n == 0 or any(len(row) != n for row in grid):
        raise ValueError("Grid must be rectangular")

    if opr_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"OPR infeasible: cost={opr_cpp['cost']} > budget={args.budget}")
    if gpr_cpp["cost"] > args.budget + 1e-9:
        raise RuntimeError(f"GPR infeasible: cost={gpr_cpp['cost']} > budget={args.budget}")
    if gpr_cpp["reward"] > opr_cpp["reward"] + 1e-9:
        raise RuntimeError(
            "Ordering violated: GPR reward exceeds OPR reward "
            f"(gpr={gpr_cpp['reward']}, opr={opr_cpp['reward']}, budget={args.budget})"
        )

    fig, axes = plt.subplots(1, 2, figsize=(max(10, n * 1.6), max(4.8, m * 0.8)), constrained_layout=True)
    fig.suptitle(
        f"Single-instance OPR vs GPR (m={m}, n={n} internal cols, graph cols={n + 2}, B={args.budget})\n"
        "Note: corridor columns are implicit and have zero profit",
        fontsize=12,
    )

    draw_instance_with_path(
        axes[0], grid, opr_cpp["cycle"], "OPR", float(opr_cpp["reward"]), float(opr_cpp["cost"]),
        None, ""
    )
    draw_instance_with_path(
        axes[1], grid, gpr_cpp["cycle"], "GPR", float(gpr_cpp["reward"]), float(gpr_cpp["cost"]),
        None, ""
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
