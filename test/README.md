`test/compare_ofr.py` compares the legacy Python OFR implementation (`_old/algorithms.py::opt_full_row`)
against the current C++ port (`src/core/algorithms.cpp`, algorithm `ofr`) on the same CSV instance.

Convention:
- the CSV stores only the `n` internal profit columns
- the two aisle-graph corridor columns (`c_0` and `c_{n+1}`) are implicit and have profit `0`

Run from repo root:

```bash
python3 test/compare_ofr.py
```

Notes:
- Requires the Release binary in `cmake-build-release/aisle-graphs`.
- The script writes a temporary shared legacy input file in `_old/input/input-990001.csv`.
