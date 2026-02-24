#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"


def compute_cycle_cost(cycle):
    if not cycle:
        return 0
    total = 0
    for (r1, c1), (r2, c2) in zip(cycle, cycle[1:]):
        total += abs(r2 - r1) + abs(c2 - c1)
    return total


def run_cpp_alg(alg_name: str, budget: int, num_rows: int, num_internal_cols: int, seed: int):
    cmd = [
        str(CPP_BIN),
        "-algorithm", alg_name,
        "-budget", str(budget),
        "-rows", str(num_rows),
        "-cols", str(num_internal_cols),
        "-seed", str(seed),
        "-log", "0",
    ]
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, capture_output=True, text=True)

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
        raise RuntimeError(
            f"Cannot parse C++ output for {alg_name}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return out


def parse_args():
    p = argparse.ArgumentParser(description="C++ APR random-instance feasibility checker.")
    p.add_argument("--rows", type=int, default=10, help="Number of rows m")
    p.add_argument("--cols", type=int, default=10, help="Number of internal columns n")
    p.add_argument("--instances", type=int, default=33, help="Number of random instances")
    p.add_argument("--seed-base", type=int, default=6000, help="Base seed (seed_base + instance_idx)")
    return p.parse_args()


def main() -> int:
    if not CPP_BIN.exists():
        raise RuntimeError(f"C++ binary not found: {CPP_BIN}")

    args = parse_args()
    num_rows = args.rows
    num_internal_cols = args.cols
    num_instances = args.instances
    seed_base = args.seed_base

    b_max = (num_internal_cols + 1) * num_rows + 2 * (num_rows - 1)
    budgets = list(range(0, b_max + 1, 2))

    print(
        f"C++ APR feasibility check on {num_instances} random instances "
        f"({num_rows}x{num_internal_cols} internal, graph cols={num_internal_cols + 2})"
    )
    print(f"Budget sweep (even only): 0..{b_max} step 2 ({len(budgets)} values)")

    total_points = 0
    violations = 0

    for inst_idx in range(num_instances):
        seed = seed_base + inst_idx
        print(f"Instance {inst_idx + 1}/{num_instances} (seed={seed})")

        for budget in budgets:
            apr = run_cpp_alg("apr", budget, num_rows, num_internal_cols, seed)
            total_points += 1

            bad = []
            if apr["cost"] > budget:
                bad.append(f"apr cost > budget ({apr['cost']} > {budget})")

            apr_path_cost = compute_cycle_cost(apr["cycle"])
            if apr_path_cost != apr["cost"]:
                bad.append(f"apr cost != path cost ({apr['cost']} vs {apr_path_cost})")

            if bad:
                violations += 1
                print(f"[FAIL] inst={inst_idx + 1} seed={seed} budget={budget}: " + "; ".join(bad))
                print(f"  apr={apr}")

    if violations:
        print(f"\nFAILED: {violations} violations over {total_points} instance-budget points.")
        return 1

    print(f"\nPASSED: {total_points} instance-budget points, 0 violations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
