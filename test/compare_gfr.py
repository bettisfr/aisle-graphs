#!/usr/bin/env python3

import csv
import math
import random
import subprocess
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
OLD_DIR = REPO_ROOT / "_old"
OLD_INPUT_DIR = OLD_DIR / "input"
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"

TEST_INPUT_ID = 990002
TEST_INPUT_FILE = OLD_INPUT_DIR / f"input-{TEST_INPUT_ID}.csv"


def write_grid_csv(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        for row in grid:
            writer.writerow(row)


def run_cpp_gfr(input_csv_path: Path, budget: int):
    if not CPP_BIN.exists():
        raise RuntimeError(
            f"C++ binary not found at {CPP_BIN}. Build it first (CLion Release / cmake-build-release)."
        )

    cmd = [
        str(CPP_BIN),
        "-algorithm", "gfr",
        "-budget", str(budget),
        "-input_csv", str(input_csv_path),
        "-log", "0",
    ]
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, capture_output=True, text=True)

    reward = None
    cost = None
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            reward = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            cost = float(line.split("=", 1)[1].strip())

    if reward is None or cost is None:
        raise RuntimeError(f"Unable to parse C++ output.\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

    return {"reward": reward, "cost": cost, "stdout": proc.stdout}


def greedy_full_row_reference(grid, budget: int):
    # Faithful to the MATLAB routine output path: only the greedy row-selection phase
    # determines returned reward/cost (the later real_tour reconstruction is not returned).
    reward_map = np.array(grid, dtype=float)
    m, n = reward_map.shape

    reward_out = 0.0
    cost_out = 0
    if budget <= 0 or m == 0 or n == 0:
        return {"reward": reward_out, "cost": cost_out}

    beginning_row = math.ceil(1 / n)  # MATLAB indexing semantics -> 1
    ending_row = math.ceil(1 / n)     # -> 1
    compiled_row_cost = n

    compiled_rewards = []
    for i in range(m):
        compiled_rewards.append((float(np.sum(reward_map[i, :])), i + 1))  # 1-based row id

    compiled_rewards.sort(key=lambda x: (-x[0], x[1]))

    current_row = beginning_row
    current_side = 1
    total_cost = 0
    total_reward = 0.0

    for row_reward, row_id in compiled_rewards:
        if current_side == 1:
            loop_cost = abs(current_row - row_id) + 2 * compiled_row_cost + abs(row_id - ending_row)
        else:
            loop_cost = abs(current_row - row_id) + compiled_row_cost + abs(row_id - ending_row)

        if loop_cost <= (budget - total_cost):
            total_cost += abs(current_row - row_id) + compiled_row_cost
            total_reward += row_reward
            current_row = row_id
            current_side = 2 if current_side == 1 else 1

    if current_side == 2:
        total_cost += abs(current_row - ending_row) + compiled_row_cost
    else:
        total_cost += abs(current_row - ending_row)

    return {"reward": float(total_reward), "cost": float(total_cost)}


def almost_equal(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def main() -> int:
    rng = random.Random(20260222)

    # Reward matrices store only the n internal profit columns.
    # Corridor columns c_0 and c_{n+1} are implicit and always have profit 0.
    fixed_cases = [
        [
            [5, 1, 2, 8, 0, 3],
            [1, 7, 4, 0, 6, 1],
            [9, 2, 1, 1, 2, 4],
            [3, 3, 8, 2, 5, 0],
            [2, 6, 0, 7, 1, 9],
        ]
    ]

    random_cases = []
    shapes = [
        (4, 6), (5, 6), (6, 8), (8, 10),
        (10, 12), (12, 16), (15, 20), (20, 24),
    ]
    for case_idx in range(48):
        m, inner_cols = shapes[case_idx % len(shapes)]
        grid = []
        for _ in range(m):
            grid.append([rng.randint(0, 10) for _ in range(inner_cols)])
        random_cases.append(grid)

    all_cases = fixed_cases + random_cases

    print(f"Testing GFR on shared instance file: {TEST_INPUT_FILE}")
    print(f"Instances: {len(all_cases)}")

    failures = 0
    comparisons = 0

    for case_no, grid in enumerate(all_cases, start=1):
        write_grid_csv(TEST_INPUT_FILE, grid)
        m = len(grid)
        inner_cols = len(grid[0])

        max_budget = inner_cols * m + 3 * m
        budgets = sorted(set(int(round(x)) for x in np.linspace(0, max_budget, num=45)))

        print(
            f"Instance {case_no}/{len(all_cases)}: "
            f"m={m}, n={inner_cols} internal cols (graph cols={inner_cols + 2}), budgets={len(budgets)}"
        )

        for budget in budgets:
            py_out = greedy_full_row_reference(grid, budget)
            cpp_out = run_cpp_gfr(TEST_INPUT_FILE, budget)
            comparisons += 1

            reward_ok = almost_equal(py_out["reward"], cpp_out["reward"])
            cost_ok = almost_equal(py_out["cost"], cpp_out["cost"])

            if not (reward_ok and cost_ok):
                failures += 1
                print(
                    f"[FAIL] case={case_no} budget={budget} "
                    f"py(reward={py_out['reward']}, cost={py_out['cost']}) "
                    f"cpp(reward={cpp_out['reward']}, cost={cpp_out['cost']})"
                )

    if failures:
        print(f"\nGFR comparison failed: {failures} mismatches out of {comparisons} comparisons.")
        return 1

    print(f"\nGFR comparison passed: {comparisons} comparisons, 0 mismatches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
