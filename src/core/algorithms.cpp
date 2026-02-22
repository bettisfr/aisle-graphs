#include "algorithms.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;

namespace {

double compute_cost_from_cycle(const vector<pair<int, int>> &cycle) {
    if (cycle.empty()) {
        return 0.0;
    }

    long long total = 0;
    for (size_t i = 0; i + 1 < cycle.size(); ++i) {
        const auto [r1, c1] = cycle[i];
        const auto [r2, c2] = cycle[i + 1];
        total += llabs(static_cast<long long>(r2) - static_cast<long long>(r1));
        total += llabs(static_cast<long long>(c2) - static_cast<long long>(c1));
    }
    return static_cast<double>(total);
}

vector<pair<int, int>> build_full_row_cycle_from_selection(const vector<int> &selected_rows_0based,
                                                           const int n_internal) {
    vector<pair<int, int>> cycle;
    if (n_internal <= 0) {
        return cycle;
    }

    const int left_corridor = 0;
    const int right_corridor = n_internal + 1;
    int side_col = left_corridor;

    cycle.push_back({0, side_col});
    for (const int row : selected_rows_0based) {
        cycle.push_back({row, side_col});
        side_col = (side_col == left_corridor) ? right_corridor : left_corridor;
        cycle.push_back({row, side_col});
    }
    cycle.push_back({0, side_col});
    cycle.push_back({0, left_corridor});

    vector<pair<int, int>> cleaned;
    cleaned.reserve(cycle.size());
    for (const auto &p : cycle) {
        if (cleaned.empty() || cleaned.back() != p) {
            cleaned.push_back(p);
        }
    }
    return cleaned;
}

} // namespace

algorithms::algorithms(const deployment &dep) : dep(dep) {
    algorithm_functions = {
        &algorithms::opt_full_row,                          // 0 -> ofr
        &algorithms::opt_partial_row_single_column,         // 1 -> oprsc
        &algorithms::opt_partial_row_single_column_novertical, // 2 -> oprsc-nv
        &algorithms::heuristic_partial_row,                 // 3 -> hprgc
        &algorithms::greedy_full_row,                       // 4 -> gfr
        &algorithms::greedy_partial_row_single_column,      // 5 -> gprsc
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
    out.cycle = build_full_row_cycle_from_selection(out.traversed_rows, n);
    out.cost = compute_cost_from_cycle(out.cycle);
    return out;
}

solution algorithms::opt_partial_row_single_column(const int budget) const {
    return make_todo_solution("oprsc", budget,
        "Port Python _old/algorithms.py::opt_partial_row_single_column (or _old/main.cpp oprsc v1)");
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
            build_full_row_cycle_from_selection(candidate_rows_0based, n);
        const double candidate_cost = compute_cost_from_cycle(candidate_cycle);

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

    out.cycle = build_full_row_cycle_from_selection(out.traversed_rows, n);
    out.cost = compute_cost_from_cycle(out.cycle);

    return out;
}

solution algorithms::greedy_partial_row_single_column(const int budget) const {
    return make_todo_solution("gprsc", budget,
        "Port greedy_partial_row_single_column baseline from Python");
}
