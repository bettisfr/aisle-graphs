#include "algorithms.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../util/util.h"

using namespace std;

algorithms::algorithms(const deployment &dep) : dep(dep) {
    algorithm_functions = {
        &algorithms::opt_full_row,                  // 0 -> ofr
        &algorithms::greedy_full_row,               // 1 -> gfr

        &algorithms::heuristic_partial_row,         // 2 -> hprgc

        &algorithms::opt_partial_row_single_column, // 3 -> oprsc
        &algorithms::greedy_partial_row_single_column, // 4 -> gprsc
    };
}

solution algorithms::run_experiment(const int algorithm, const int budget) {
    solution out;

    if (algorithm >= 0 && algorithm < static_cast<int>(algorithm_functions.size())) {
        out = algorithm_functions[algorithm](*this, budget);
    } else {
        throw invalid_argument("Invalid algorithm index: " + to_string(algorithm));
    }

    if (out.cost > static_cast<double>(budget) + 1e-9) {
        throw runtime_error(
            "Infeasible solution: algorithm cost=" + to_string(out.cost) +
            " exceeds budget=" + to_string(budget)
        );
    }

    // Validate full-row solutions structurally (path/layout/cost consistency).
    if (algorithm == 0 || algorithm == 1) { // ofr, gfr
        string reason;
        if (!util::is_full_row_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid full-row solution: " + reason);
        }
    }
    if (algorithm == 3 || algorithm == 4) { // oprsc, gprsc
        string reason;
        if (!util::is_partial_row_single_column_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid partial-row single-column solution: " + reason);
        }
    }

    return out;
}

solution algorithms::make_todo_solution(const string &algorithm_key, const int budget, const string &note) const {
    solution out;
    out.algorithm_key = algorithm_key;
    out.reward = 0.0;
    out.cost = 0.0;
    out.notes.push_back("TODO: " + note);
    out.notes.push_back("Budget received=" + to_string(budget));
    out.notes.push_back("Instance dims=" + to_string(dep.rows()) + "x" + to_string(dep.cols()));
    return out;
}

solution algorithms::opt_full_row(const int budget) const {
    solution out;
    out.algorithm_key = "ofr";

    const RewardGrid &rewards = dep.rewards();
    const int m = dep.rows();
    const int n = dep.cols();

    if (budget <= 0 || m <= 0 || n <= 0) {
        return out;
    }

    // OPTFRI-style implementation under the unified cost metric:
    // choose an optimal subset of full rows; cost depends only on max selected row and number of row crossings.
    // Each full-row crossing on the full graph costs (n+1) edges (from corridor to corridor).
    const int row_cross_cost = n + 1;

    vector<double> row_rewards(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            row_rewards[i] += rewards[i][j];
        }
    }

    double best_reward = 0.0;
    double best_cost = 0.0;
    vector<int> best_rows;

    for (int max_row = 0; max_row < m; ++max_row) {
        // Candidate subsets must include max_row and can only contain rows in [0..max_row].
        vector<int> prefix_ids;
        prefix_ids.reserve(max_row);
        for (int i = 0; i < max_row; ++i) {
            prefix_ids.push_back(i);
        }
        sort(prefix_ids.begin(), prefix_ids.end(),
             [&row_rewards](const int a, const int b) {
                 if (row_rewards[a] == row_rewards[b]) {
                     return a < b;
                 }
                 return row_rewards[a] > row_rewards[b];
             });

        // Try all possible numbers of selected rows k >= 1 including max_row.
        for (int k = 1; k <= max_row + 1; ++k) {
            // Unified full-graph cost:
            // vertical down/up to deepest row + k row crossings + optional final top-row crossing if k is odd.
            const int crossings = k + (k % 2);
            const double candidate_cost =
                static_cast<double>(2 * max_row + crossings * row_cross_cost);

            if (candidate_cost > static_cast<double>(budget)) {
                continue;
            }

            vector<int> candidate_rows;
            candidate_rows.reserve(k);

            // Select k-1 best rows among [0..max_row-1], then include max_row.
            for (int t = 0; t < k - 1 && t < static_cast<int>(prefix_ids.size()); ++t) {
                candidate_rows.push_back(prefix_ids[t]);
            }
            candidate_rows.push_back(max_row);

            sort(candidate_rows.begin(), candidate_rows.end());

            double candidate_reward = 0.0;
            for (const int r : candidate_rows) {
                candidate_reward += row_rewards[r];
            }

            if (candidate_reward > best_reward) {
                best_reward = candidate_reward;
                best_cost = candidate_cost;
                best_rows = candidate_rows;
            }
        }
    }

    if (best_rows.empty()) {
        return out;
    }

    out.reward = best_reward;
    out.traversed_rows = best_rows;
    out.cycle = util::build_full_row_cycle_from_selection(out.traversed_rows, n);
    out.cost = util::compute_cycle_cost(out.cycle);
    return out;
}

solution algorithms::opt_partial_row_single_column(const int budget) const {
    solution out;
    out.algorithm_key = "oprsc";

    const RewardGrid &rewards = dep.rewards();
    const int m = dep.rows();
    const int n = dep.cols(); // internal profitable columns

    if (budget <= 0 || m <= 0 || n <= 0) {
        return out;
    }

    // Prefix rewards by number of visited internal columns:
    // pref[i][c] = reward collected on row i by visiting columns 1..c (full-graph coords), with c in [0..n].
    vector<vector<double>> pref(m, vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int c = 1; c <= n; ++c) {
            pref[i][c] = pref[i][c - 1] + rewards[i][c - 1];
        }
    }

    // Reduced budget units: every movement contribution is even under this policy
    // (vertical on left corridor + out-and-back horizontal), so we work on B/2.
    const int B2 = budget / 2;
    if (B2 <= 0) {
        out.cycle = {{0, 0}};
        return out;
    }

    const double NEG_INF = -1e100;
    vector<vector<double>> R(m, vector<double>(B2 + 1, NEG_INF));
    vector<vector<int>> S(m, vector<int>(B2 + 1, 0)); // chosen prefix length c in [0..n]

    // First row: no vertical cost, only horizontal out-and-back cost 2*c.
    for (int b = 0; b <= B2; ++b) {
        const int best_c = min(n, b);
        R[0][b] = pref[0][best_c];
        S[0][b] = best_c;
    }

    // DP on rows. Transition consumes c + 1 reduced units:
    // +1 accounts for moving the "active frontier" one row deeper.
    for (int i = 1; i < m; ++i) {
        for (int b = 0; b <= B2; ++b) {
            if (b < i) {
                R[i][b] = NEG_INF;
                S[i][b] = 0;
                continue;
            }

            if (b == i) {
                // Reach row i and come back with zero horizontal exploration on all rows.
                R[i][b] = 0.0;
                S[i][b] = 0;
                continue;
            }

            double best_val = NEG_INF;
            int best_c = 0;
            for (int c = 0; c <= n; ++c) {
                const int idx = b - c - 1;
                if (idx < 0) {
                    continue;
                }
                if (R[i - 1][idx] <= NEG_INF / 2) {
                    continue;
                }
                const double cand = R[i - 1][idx] + pref[i][c];
                if (cand > best_val) {
                    best_val = cand;
                    best_c = c;
                }
            }

            R[i][b] = best_val;
            S[i][b] = best_c;
        }
    }

    // Pick best deepest visited row at full reduced budget.
    int best_max_row = -1;
    double best_reward = 0.0;
    for (int i = 0; i < m; ++i) {
        if (R[i][B2] > best_reward) {
            best_reward = R[i][B2];
            best_max_row = i;
        }
    }

    if (best_max_row < 0 || best_reward <= 0.0) {
        out.cycle = {{0, 0}};
        return out;
    }

    // Backtrack selected prefix lengths c_i for rows 0..best_max_row.
    vector<int> c_by_row(best_max_row + 1, 0);
    int b = B2;
    for (int i = best_max_row; i >= 0; --i) {
        const int c = S[i][b];
        c_by_row[i] = c;
        if (i > 0) {
            b -= (c + 1);
            if (b < 0) {
                throw runtime_error("oprsc backtracking produced negative reduced budget");
            }
        }
    }

    // Build path on full graph coordinates (corridors included): only left corridor verticals.
    vector<pair<int, int>> cycle;
    cycle.reserve(3 * (best_max_row + 1) + 2);
    cycle.push_back({0, 0});
    for (int i = 0; i <= best_max_row; ++i) {
        if (cycle.back() != make_pair(i, 0)) {
            cycle.push_back({i, 0});
        }
        if (c_by_row[i] > 0) {
            cycle.push_back({i, c_by_row[i]});
            cycle.push_back({i, 0});
            out.traversed_rows.push_back(i);
        }
        out.selected_columns_per_row.push_back(c_by_row[i]);
    }
    if (cycle.back() != make_pair(0, 0)) {
        cycle.push_back({0, 0});
    }

    out.reward = best_reward;
    out.cycle = move(cycle);
    out.cost = util::compute_cycle_cost(out.cycle);
    return out;
}

solution algorithms::opt_partial_row_single_column_novertical(const int budget) const {
    return make_todo_solution("oprsc-nv", budget,
        "Port no-vertical DP variant from _old/algorithms.py::opt_partial_row_single_column_novertical / _old/main.cpp v2");
}

solution algorithms::heuristic_partial_row(const int budget) const {
    return make_todo_solution("hprgc", budget,
        "Port heuristic_partial_row + heuristic_1..5 after OFR/OPRSC are available");
}

solution algorithms::greedy_full_row(const int budget) const {
    solution out;
    out.algorithm_key = "gfr";

    const RewardGrid &reward_map = dep.rewards();
    const int m = dep.rows();
    const int n = dep.cols(); // internal columns (same matrix convention as legacy code)

    if (budget <= 0 || m <= 0 || n <= 0) {
        return out;
    }

    vector<pair<double, int>> rows; // {row_reward, row_id_1based}
    rows.reserve(m);
    for (int i = 0; i < m; ++i) {
        double row_reward = 0.0;
        for (int j = 0; j < n; ++j) {
            row_reward += reward_map[i][j];
        }
        rows.emplace_back(row_reward, i + 1);
    }

    sort(rows.begin(), rows.end(),
         [](const pair<double, int> &a, const pair<double, int> &b) {
             if (a.first == b.first) {
                 return a.second < b.second;
             }
             return a.first > b.first;
         });

    double total_reward = 0.0;
    vector<int> selected_rows_1based;

    for (const int row_id : selected_rows_1based) {
        (void) row_id;
    }
    selected_rows_1based.clear();

    for (const auto &[row_reward, row_id] : rows) {
        vector<int> candidate_rows_0based;
        candidate_rows_0based.reserve(selected_rows_1based.size() + 1);
        for (const int rid : selected_rows_1based) {
            candidate_rows_0based.push_back(rid - 1);
        }
        candidate_rows_0based.push_back(row_id - 1);

        const vector<pair<int, int>> candidate_cycle =
            util::build_full_row_cycle_from_selection(candidate_rows_0based, n);
        const double candidate_cost = util::compute_cycle_cost(candidate_cycle);

        if (candidate_cost <= static_cast<double>(budget)) {
            selected_rows_1based.push_back(row_id);
            total_reward += row_reward;
        }
    }

    out.reward = total_reward;
    out.traversed_rows.reserve(selected_rows_1based.size());
    for (const int row_id : selected_rows_1based) {
        out.traversed_rows.push_back(row_id - 1);
    }

    out.cycle = util::build_full_row_cycle_from_selection(out.traversed_rows, n);
    out.cost = util::compute_cycle_cost(out.cycle);

    return out;
}

solution algorithms::greedy_partial_row_single_column(const int budget) const {
    solution out;
    out.algorithm_key = "gprsc";

    const RewardGrid &rewards = dep.rewards();
    const int m = dep.rows();
    const int n = dep.cols(); // internal profitable columns

    if (budget <= 0 || m <= 0 || n <= 0) {
        return out;
    }

    // Prefix rewards on internal columns: pref[i][c] is reward from first c internal columns (c in [0..n]).
    vector<vector<double>> pref(m, vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int c = 1; c <= n; ++c) {
            pref[i][c] = pref[i][c - 1] + rewards[i][c - 1];
        }
    }

    int current_row = 0; // always on left corridor
    double total_cost = 0.0;
    double total_reward = 0.0;

    // selected_prefix[i] = number of internal columns already collected on row i (0..n)
    vector<int> selected_prefix(m, 0);
    vector<pair<int, int>> cycle;
    cycle.push_back({0, 0});

    while (true) {
        double best_score = -1.0;
        int best_row = -1;
        int best_c = 0;

        for (int i = 0; i < m; ++i) {
            for (int c = selected_prefix[i] + 1; c <= n; ++c) {
                const double add_reward = pref[i][c] - pref[i][selected_prefix[i]];
                if (add_reward <= 0.0) {
                    continue;
                }

                const double add_cost = static_cast<double>(abs(i - current_row) + 2 * c);
                const double return_home_cost = static_cast<double>(i);
                if (total_cost + add_cost + return_home_cost > static_cast<double>(budget) + 1e-9) {
                    continue;
                }

                if (add_cost <= 0.0) {
                    continue;
                }

                const double score = add_reward / add_cost;
                if (score > best_score + 1e-12) {
                    best_score = score;
                    best_row = i;
                    best_c = c;
                } else if (fabs(score - best_score) <= 1e-12 && best_row != -1) {
                    // Deterministic tie-break close to legacy numpy argmax scanning order.
                    if (i < best_row || (i == best_row && c < best_c)) {
                        best_row = i;
                        best_c = c;
                    }
                }
            }
        }

        if (best_row < 0) {
            break;
        }

        // Move vertically on the left corridor to the chosen row.
        total_cost += static_cast<double>(abs(best_row - current_row));
        if (cycle.back() != make_pair(best_row, 0)) {
            cycle.push_back({best_row, 0});
        }

        // Partial-row excursion on the full graph: to internal column 'best_c' and back to left corridor.
        cycle.push_back({best_row, best_c});
        cycle.push_back({best_row, 0});
        total_cost += static_cast<double>(2 * best_c);

        total_reward += pref[best_row][best_c] - pref[best_row][selected_prefix[best_row]];
        selected_prefix[best_row] = best_c;
        out.traversed_rows.push_back(best_row);
        current_row = best_row;
    }

    if (current_row != 0) {
        cycle.push_back({0, 0});
        total_cost += static_cast<double>(current_row);
    }

    out.reward = total_reward;
    out.selected_columns_per_row = selected_prefix;
    out.cycle = move(cycle);
    out.cost = util::compute_cycle_cost(out.cycle);

    // Keep reward/cost based on the explicit path/selection we built.
    // The path-cost recomputation should match total_cost; keep a note if not for debugging.
    if (fabs(out.cost - total_cost) > 1e-9) {
        out.notes.push_back("Warning: internal cost accumulator mismatch with cycle cost");
    }

    return out;
}
