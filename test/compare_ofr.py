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

# Large id to avoid collisions with legacy datasets.
TEST_INPUT_ID = 990001
TEST_INPUT_FILE = OLD_INPUT_DIR / f"input-{TEST_INPUT_ID}.csv"


def _install_win32com_stub() -> None:
    # algorithms.py imports win32com.client unconditionally, but OFR does not use it.
    if "win32com.client" in sys.modules:
        return
    win32com = types.ModuleType("win32com")
    client = types.ModuleType("win32com.client")
    win32com.client = client
    sys.modules["win32com"] = win32com
    sys.modules["win32com.client"] = client


def _install_legacy_util_stub() -> None:
    # algorithms.py imports `from util import *`; for OFR we only need `read_input_from_file_id`.
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


def write_grid_csv(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        for row in grid:
            writer.writerow(row)


def load_python_ofr():
    _install_win32com_stub()
    _install_legacy_util_stub()
    sys.path.insert(0, str(OLD_DIR))
    import algorithms as old_algorithms  # type: ignore
    return old_algorithms.opt_full_row


def run_cpp_ofr(input_csv_path: Path, budget: int):
    if not CPP_BIN.exists():
        raise RuntimeError(
            f"C++ binary not found at {CPP_BIN}. Build it first (CLion Release / cmake-build-release)."
        )

    cmd = [
        str(CPP_BIN),
        "-algorithm", "ofr",
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


def almost_equal(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def main() -> int:
    opt_full_row_py = load_python_ofr()

    rng = random.Random(123456)

    # Reward matrices store only the n internal profit columns.
    # Corridor columns c_0 and c_{n+1} are implicit and always have profit 0.
    # First case fixed and easy to inspect, then many random cases.
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
            row = [rng.randint(0, 10) for _ in range(inner_cols)]
            grid.append(row)
        random_cases.append(grid)

    all_cases = fixed_cases + random_cases

    print(f"Testing OFR on shared instance file: {TEST_INPUT_FILE}")
    print(f"Instances: {len(all_cases)}")

    failures = 0
    comparisons = 0

    for case_no, grid in enumerate(all_cases, start=1):
        write_grid_csv(TEST_INPUT_FILE, grid)
        m = len(grid)
        inner_cols = len(grid[0])

        max_budget = inner_cols * m + 3 * m  # same spirit as execute_batch()
        budgets = sorted(set(int(round(x)) for x in np.linspace(0, max_budget, num=35)))

        print(
            f"Instance {case_no}/{len(all_cases)}: "
            f"m={m}, n={inner_cols} internal cols (graph cols={inner_cols + 2}), budgets={len(budgets)}"
        )

        for budget in budgets:
            py_out = opt_full_row_py(TEST_INPUT_ID, budget)
            cpp_out = run_cpp_ofr(TEST_INPUT_FILE, budget)
            comparisons += 1

            reward_ok = almost_equal(float(py_out["reward"]), float(cpp_out["reward"]))
            cost_ok = almost_equal(float(py_out["cost"]), float(cpp_out["cost"]))

            if not (reward_ok and cost_ok):
                failures += 1
                print(
                    f"[FAIL] case={case_no} budget={budget} "
                    f"py(reward={py_out['reward']}, cost={py_out['cost']}) "
                    f"cpp(reward={cpp_out['reward']}, cost={cpp_out['cost']})"
                )

    if failures:
        print(f"\nOFR comparison failed: {failures} mismatches out of {comparisons} comparisons.")
        return 1

    print(f"\nOFR comparison passed: {comparisons} comparisons, 0 mismatches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
