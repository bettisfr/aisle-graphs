#!/usr/bin/env python3

import argparse
import math
import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"

ALGORITHMS = ["opr", "apr", "hpr", "gpr"]
OPR_MAX_BUDGET = 40
COLORS = {
    "opr": "#1f77b4",
    "apr": "#ff7f0e",
    "hpr": "#2ca02c",
    "gpr": "#d62728",
}


def compute_cycle_cost(cycle):
    if not cycle:
        return 0
    total = 0
    for (r1, c1), (r2, c2) in zip(cycle, cycle[1:]):
        total += abs(r2 - r1) + abs(c2 - c1)
    return total


def run_cpp_alg(alg_name: str, budget: int, rows: int, cols: int, seed: int, env: dict) -> dict:
    cmd = [
        str(CPP_BIN),
        "-algorithm", alg_name,
        "-budget", str(budget),
        "-rows", str(rows),
        "-cols", str(cols),
        "-seed", str(seed),
        "-log", "0",
    ]
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=False, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(
            f"{alg_name} failed at seed={seed}, budget={budget}, returncode={proc.returncode}. "
            f"Message: {msg}"
        )
    out = {"reward": None, "cost": None, "cycle": []}
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            out["reward"] = int(float(line.split("=", 1)[1].strip()))
        elif line.startswith("cost="):
            out["cost"] = int(float(line.split("=", 1)[1].strip()))
        elif line.startswith("cycle_list="):
            raw = line.split("=", 1)[1].strip()
            if raw:
                for token in raw.split(";"):
                    if not token:
                        continue
                    r_str, c_str = token.split(":")
                    out["cycle"].append((int(r_str), int(c_str)))

    if out["reward"] is None or out["cost"] is None:
        raise RuntimeError(f"Cannot parse output for {alg_name} at B={budget}, seed={seed}")

    if out["cost"] > budget:
        raise RuntimeError(
            f"Infeasible cost for {alg_name}: cost={out['cost']} > budget={budget} (seed={seed})"
        )

    if not out["cycle"]:
        if out["cost"] != 0 or out["reward"] != 0:
            raise RuntimeError(
                f"Empty cycle with non-zero cost/reward for {alg_name} at B={budget}, seed={seed}"
            )
        return out

    if out["cycle"][0] != (0, 0) or out["cycle"][-1] != (0, 0):
        raise RuntimeError(f"Cycle must start/end at (0,0) for {alg_name} at B={budget}, seed={seed}")

    path_cost = compute_cycle_cost(out["cycle"])
    if path_cost != out["cost"]:
        raise RuntimeError(
            f"Path-cost mismatch for {alg_name}: reported={out['cost']} path={path_cost} "
            f"(B={budget}, seed={seed})"
        )

    return out


def parse_args():
    p = argparse.ArgumentParser(
        description="Evaluate PR algorithms (OPR/APR/HPR/GPR) on random instances and plot average reward vs budget."
    )
    p.add_argument("--rows", type=int, default=6, help="Number of rows m")
    p.add_argument("--cols", type=int, default=6, help="Number of internal columns n")
    p.add_argument("--instances", type=int, default=33, help="Number of random instances")
    p.add_argument("--seed-base", type=int, default=7000, help="Base seed (seed_base + instance_idx)")
    p.add_argument("--budget-step", type=int, default=1, help="Budget step for sweep")
    p.add_argument("--min-budget", type=int, default=0, help="Minimum budget for sweep")
    p.add_argument("--max-budget", type=int, default=-1, help="Optional max budget override (<= computed Bmax)")
    p.add_argument("--out", type=str, default="test/output/evaluate_pr.pdf", help="Output PDF path")
    p.add_argument("--grb-license-file", type=str, default="", help="Optional GRB_LICENSE_FILE value")
    p.add_argument("--disable-opr", action="store_true", help="Do not execute/plot OPR")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not CPP_BIN.exists():
        raise RuntimeError(f"C++ binary not found: {CPP_BIN}")
    if args.budget_step <= 0:
        raise ValueError("--budget-step must be > 0")
    if args.min_budget < 0:
        raise ValueError("--min-budget must be >= 0")

    out_path = Path(args.out)
    if not out_path.parent.exists():
        raise RuntimeError(f"Output directory does not exist: {out_path.parent}")
    if not out_path.parent.is_dir():
        raise RuntimeError(f"Output path parent is not a directory: {out_path.parent}")

    env = os.environ.copy()
    if args.grb_license_file:
        env["GRB_LICENSE_FILE"] = args.grb_license_file
    active_algorithms = [alg for alg in ALGORITHMS if not (args.disable_opr and alg == "opr")]

    # Paper upper bound using internal columns n convention.
    b_max = (args.cols + 1) * args.rows + 2 * (args.rows - 1)
    if args.max_budget >= 0:
        b_max = min(b_max, args.max_budget)
    budgets = [b for b in range(args.min_budget, b_max + 1, args.budget_step) if b % 2 == 0]
    if not budgets:
        raise ValueError("No even budgets generated: adjust --budget-step/--max-budget")

    values = {alg: [[] for _ in budgets] for alg in active_algorithms}

    print(
        f"Evaluating PR algorithms on {args.instances} random instances "
        f"({args.rows}x{args.cols} internal, graph cols={args.cols + 2})"
    )
    print(f"Budget sweep: {budgets[0]}..{budgets[-1]} step {args.budget_step} ({len(budgets)} values)")

    for inst_idx in range(args.instances):
        seed = args.seed_base + inst_idx
        print(f"Instance {inst_idx + 1}/{args.instances} (seed={seed})")
        for bi, budget in enumerate(budgets):
            by_alg = {}
            for alg in active_algorithms:
                if alg == "opr" and budget > OPR_MAX_BUDGET:
                    continue
                by_alg[alg] = run_cpp_alg(alg, budget, args.rows, args.cols, seed, env)

            if "opr" in by_alg:
                opr_reward = by_alg["opr"]["reward"]
                # OPR is the exact internal baseline; enforce dominance only when OPR is executed.
                for alg in ["apr", "hpr", "gpr"]:
                    if by_alg[alg]["reward"] > opr_reward:
                        raise RuntimeError(
                            "Ordering violated: "
                            f"reward(opr)={opr_reward} < reward({alg})={by_alg[alg]['reward']} "
                            f"(seed={seed}, budget={budget})"
                        )

            for alg in active_algorithms:
                if alg in by_alg:
                    values[alg][bi].append(by_alg[alg]["reward"])

    means = {}
    stds = {}
    for alg in active_algorithms:
        alg_means = []
        alg_stds = []
        for arr in values[alg]:
            if not arr:
                alg_means.append(float("nan"))
                alg_stds.append(float("nan"))
                continue
            mu = sum(arr) / len(arr)
            if len(arr) <= 1:
                sigma = 0.0
            else:
                var = sum((x - mu) ** 2 for x in arr) / (len(arr) - 1)
                sigma = math.sqrt(var)
            alg_means.append(mu)
            alg_stds.append(sigma)
        means[alg] = alg_means
        stds[alg] = alg_stds

    plt.figure(figsize=(9, 5.2))
    for alg in active_algorithms:
        mean = means[alg]
        std = stds[alg]
        lower = [m - s for m, s in zip(mean, std)]
        upper = [m + s for m, s in zip(mean, std)]
        plt.plot(budgets, mean, label=alg.upper(), linewidth=2.0, color=COLORS[alg])
        plt.fill_between(budgets, lower, upper, color=COLORS[alg], alpha=0.15)

    plt.title(
        f"PR Algorithms on Random {args.rows}x{args.cols} Instances "
        f"({args.instances} seeds)"
    )
    plt.xlabel("Budget")
    plt.ylabel("Reward (mean ± std)")
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    print(f"Saved figure to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
