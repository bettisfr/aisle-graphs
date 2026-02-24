#!/usr/bin/env python3

import csv
import random
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"
TMP_INPUT = REPO_ROOT / "test" / "output" / "tmp_check_cpp_ofr_gfr_random_10x10.csv"


def write_grid_csv(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        csv.writer(f).writerows(grid)


def run_cpp_alg(alg_name: str, input_csv_path: Path, budget: int):
    cmd = [
        str(CPP_BIN),
        "-algorithm", alg_name,
        "-budget", str(budget),
        "-input_csv", str(input_csv_path),
        "-log", "0",
    ]
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, capture_output=True, text=True)

    out = {
        "reward": None,
        "cost": None,
        "full_row_feasible": None,
        "full_row_reason": "",
    }
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            out["reward"] = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            out["cost"] = float(line.split("=", 1)[1].strip())
        elif line.startswith("full_row_feasible="):
            out["full_row_feasible"] = int(line.split("=", 1)[1].strip())
        elif line.startswith("full_row_reason="):
            out["full_row_reason"] = line.split("=", 1)[1].strip()

    if out["reward"] is None or out["cost"] is None:
        raise RuntimeError(
            f"Cannot parse C++ output for {alg_name}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return out


def build_random_grid(rows: int, cols: int, seed: int):
    rng = random.Random(seed)
    return [[rng.randint(0, 10) for _ in range(cols)] for _ in range(rows)]


def main() -> int:
    if not CPP_BIN.exists():
        raise RuntimeError(f"C++ binary not found: {CPP_BIN}")

    rows = 10
    internal_cols = 10
    num_instances = 33

    # Paper upper bound for exhaustive visit feasibility:
    # B_max = (n + 1)m + 2(m - 1), where n is the number of internal columns.
    b_max = (internal_cols + 1) * rows + 2 * (rows - 1)
    budgets = list(range(0, b_max + 1))

    print(
        f"C++ OFR/GFR check on {num_instances} random instances "
        f"({rows}x{internal_cols} internal, graph cols={internal_cols + 2})"
    )
    print(f"Budget sweep: 0..{b_max} ({len(budgets)} values)")

    total_points = 0
    violations = 0

    for inst_idx in range(num_instances):
        seed = 1000 + inst_idx
        grid = build_random_grid(rows, internal_cols, seed)
        write_grid_csv(TMP_INPUT, grid)

        print(f"Instance {inst_idx + 1}/{num_instances} (seed={seed})")

        for budget in budgets:
            ofr = run_cpp_alg("ofr", TMP_INPUT, budget)
            gfr = run_cpp_alg("gfr", TMP_INPUT, budget)
            total_points += 1

            bad = []
            if ofr["cost"] > budget + 1e-9:
                bad.append(f"ofr cost > budget ({ofr['cost']} > {budget})")
            if gfr["cost"] > budget + 1e-9:
                bad.append(f"gfr cost > budget ({gfr['cost']} > {budget})")
            if ofr["full_row_feasible"] != 1:
                bad.append(f"ofr validator KO ({ofr['full_row_reason']})")
            if gfr["full_row_feasible"] != 1:
                bad.append(f"gfr validator KO ({gfr['full_row_reason']})")
            if gfr["reward"] > ofr["reward"] + 1e-9:
                bad.append(f"ordering violated gfr > ofr ({gfr['reward']} > {ofr['reward']})")

            if bad:
                violations += 1
                print(
                    f"[FAIL] inst={inst_idx + 1} seed={seed} budget={budget}: " + "; ".join(bad)
                )
                print(f"  ofr={ofr}")
                print(f"  gfr={gfr}")

    if violations:
        print(f"\nFAILED: {violations} violations over {total_points} instance-budget points.")
        return 1

    print(f"\nPASSED: {total_points} instance-budget points, 0 violations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
