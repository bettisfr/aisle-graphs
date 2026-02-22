#!/usr/bin/env python3

import csv
import math
import random
import subprocess
import sys
import types
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
OLD_DIR = REPO_ROOT / "_old"
OLD_INPUT_DIR = OLD_DIR / "input"
CPP_BIN = REPO_ROOT / "cmake-build-release" / "aisle-graphs"

TEST_INPUT_ID = 990003
TEST_INPUT_FILE = OLD_INPUT_DIR / f"input-{TEST_INPUT_ID}.csv"


def _install_win32com_stub() -> None:
    if "win32com.client" in sys.modules:
        return
    win32com = types.ModuleType("win32com")
    client = types.ModuleType("win32com.client")
    win32com.client = client
    sys.modules["win32com"] = win32com
    sys.modules["win32com.client"] = client


def _install_legacy_util_stub() -> None:
    if "util" in sys.modules:
        return

    util_mod = types.ModuleType("util")

    def read_input_from_file_id(input_id):
        path = OLD_INPUT_DIR / f"input-{int(input_id)}.csv"
        with path.open("r", newline="") as f:
            reader = csv.reader(f, delimiter=",")
            x = list(reader)
        return np.array(x).astype("float")

    util_mod.read_input_from_file_id = read_input_from_file_id
    sys.modules["util"] = util_mod


def load_legacy_ofr():
    _install_win32com_stub()
    _install_legacy_util_stub()
    if str(OLD_DIR) not in sys.path:
        sys.path.insert(0, str(OLD_DIR))
    import algorithms as old_algorithms  # type: ignore
    return old_algorithms.opt_full_row


def greedy_full_row_reference(grid, budget: int):
    reward_map = np.array(grid, dtype=float)
    m, n = reward_map.shape

    if budget <= 0 or m == 0 or n == 0:
        return {"reward": 0.0, "cost": 0.0}

    beginning_row = math.ceil(1 / n)
    ending_row = math.ceil(1 / n)
    compiled_row_cost = n

    compiled_rewards = []
    for i in range(m):
        compiled_rewards.append((float(np.sum(reward_map[i, :])), i + 1))
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


def write_grid_csv(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(grid)


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
    for line in proc.stdout.splitlines():
        if line.startswith("reward="):
            reward = float(line.split("=", 1)[1].strip())
        elif line.startswith("cost="):
            cost = float(line.split("=", 1)[1].strip())

    if reward is None or cost is None:
        raise RuntimeError(f"Cannot parse output for {alg_name}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

    return {"reward": reward, "cost": cost}


def almost_equal(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def main() -> int:
    legacy_ofr = load_legacy_ofr()
    rng = random.Random(424242)

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
    shapes = [(4, 6), (5, 6), (6, 8), (8, 10), (10, 12), (12, 16)]
    for case_idx in range(36):
        m, inner_cols = shapes[case_idx % len(shapes)]
        random_cases.append([[rng.randint(0, 10) for _ in range(inner_cols)] for _ in range(m)])

    all_cases = fixed_cases + random_cases

    print(f"Testing OFR/GFR consistency on shared instance file: {TEST_INPUT_FILE}")
    print(f"Instances: {len(all_cases)}")

    failures = 0
    comparisons = 0

    for case_no, grid in enumerate(all_cases, start=1):
        write_grid_csv(TEST_INPUT_FILE, grid)
        m = len(grid)
        inner_cols = len(grid[0])
        max_budget = inner_cols * m + 3 * m
        budgets = sorted(set(int(round(x)) for x in np.linspace(0, max_budget, num=31)))

        print(
            f"Instance {case_no}/{len(all_cases)}: "
            f"m={m}, n={inner_cols} internal cols (graph cols={inner_cols + 2}), budgets={len(budgets)}"
        )

        for budget in budgets:
            old_ofr = legacy_ofr(TEST_INPUT_ID, budget)
            old_gfr = greedy_full_row_reference(grid, budget)
            new_ofr = run_cpp_alg("ofr", TEST_INPUT_FILE, budget)
            new_gfr = run_cpp_alg("gfr", TEST_INPUT_FILE, budget)
            comparisons += 1

            checks = [
                ("ofr reward", almost_equal(float(old_ofr["reward"]), new_ofr["reward"])),
                ("ofr cost", almost_equal(float(old_ofr["cost"]), new_ofr["cost"])),
                ("gfr reward", almost_equal(float(old_gfr["reward"]), new_gfr["reward"])),
                ("gfr cost", almost_equal(float(old_gfr["cost"]), new_gfr["cost"])),
                ("old ofr feasible", float(old_ofr["cost"]) <= float(budget) + 1e-9),
                ("new ofr feasible", float(new_ofr["cost"]) <= float(budget) + 1e-9),
                ("old gfr feasible", float(old_gfr["cost"]) <= float(budget) + 1e-9),
                ("new gfr feasible", float(new_gfr["cost"]) <= float(budget) + 1e-9),
                ("old gfr<=ofr", float(old_gfr["reward"]) <= float(old_ofr["reward"]) + 1e-9),
                ("new gfr<=ofr", float(new_gfr["reward"]) <= float(new_ofr["reward"]) + 1e-9),
            ]

            bad = [name for name, ok in checks if not ok]
            if bad:
                failures += 1
                print(
                    f"[FAIL] case={case_no} budget={budget} failed={bad}\n"
                    f"  old_ofr={old_ofr}\n"
                    f"  new_ofr={new_ofr}\n"
                    f"  old_gfr={old_gfr}\n"
                    f"  new_gfr={new_gfr}"
                )

    if failures:
        print(f"\nOFR/GFR consistency failed: {failures} mismatches out of {comparisons} budget-points.")
        return 1

    print(f"\nOFR/GFR consistency passed: {comparisons} budget-points, 0 mismatches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
