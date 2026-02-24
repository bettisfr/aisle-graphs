#include "algorithms.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../util/util.h"

using namespace std;

solution algorithms::solve_oprsc_left_on_grid(const RewardGrid &rewards, const int budget, const string &algorithm_key,
                                              const bool disable_first_row_partial) const {
    solution out;
    out.algorithm_key = algorithm_key;

    const int num_rows = static_cast<int>(rewards.size());
    const int num_internal_cols = (num_rows > 0) ? static_cast<int>(rewards[0].size()) : 0;

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        return out;
    }

    vector<vector<int>> pref(num_rows, vector<int>(num_internal_cols + 1, 0));
    for (int i = 0; i < num_rows; ++i) {
        for (int c = 1; c <= num_internal_cols; ++c) {
            pref[i][c] = pref[i][c - 1] + rewards[i][c - 1];
        }
    }

    const int B2 = budget / 2;
    if (B2 <= 0) {
        out.cycle = {{0, 0}};
        return out;
    }

    const int NEG_INF = -1000000000;
    vector<vector<int>> R(num_rows, vector<int>(B2 + 1, NEG_INF));
    vector<vector<int>> S(num_rows, vector<int>(B2 + 1, 0));

    for (int b = 0; b <= B2; ++b) {
        if (disable_first_row_partial) {
            R[0][b] = 0;
            S[0][b] = 0;
        } else {
            const int best_c = min(num_internal_cols, b);
            R[0][b] = pref[0][best_c];
            S[0][b] = best_c;
        }
    }

    for (int i = 1; i < num_rows; ++i) {
        for (int b = 0; b <= B2; ++b) {
            if (b < i) {
                R[i][b] = NEG_INF;
                S[i][b] = 0;
                continue;
            }
            if (b == i) {
                R[i][b] = 0;
                S[i][b] = 0;
                continue;
            }

            int best_val = NEG_INF;
            int best_c = 0;
            for (int c = 0; c <= num_internal_cols; ++c) {
                const int idx = b - c - 1;
                if (idx < 0) {
                    continue;
                }
                if (R[i - 1][idx] <= NEG_INF / 2) {
                    continue;
                }
                const int cand = R[i - 1][idx] + pref[i][c];
                if (cand > best_val) {
                    best_val = cand;
                    best_c = c;
                } else if (cand == best_val && c < best_c) {
                    // Prefer shorter partial traversal when reward is identical.
                    best_c = c;
                }
            }
            R[i][b] = best_val;
            S[i][b] = best_c;
        }
    }

    int best_max_row = -1;
    int best_reward = 0;
    for (int i = 0; i < num_rows; ++i) {
        if (R[i][B2] > best_reward) {
            best_reward = R[i][B2];
            best_max_row = i;
        }
    }

    if (best_max_row < 0 || best_reward <= 0) {
        out.cycle = {{0, 0}};
        return out;
    }

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

    // Trim useless zero-reward tails from each selected prefix so the reconstructed path is tight.
    for (int i = 0; i <= best_max_row; ++i) {
        while (c_by_row[i] > 0 && fabs(rewards[i][c_by_row[i] - 1]) <= 1e-12) {
            c_by_row[i]--;
        }
    }

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

algorithms::algorithms(const deployment &dep) : dep(dep) {
    algorithm_functions = {
        &algorithms::opt_full_row,                  // 0 -> ofr
        &algorithms::greedy_full_row,               // 1 -> gfr

        &algorithms::heuristic_partial_row,         // 2 -> hpr
        &algorithms::apx_partial_row,               // 3 -> apr

        &algorithms::opt_partial_row_single_column, // 4 -> oprsc
        &algorithms::greedy_partial_row_single_column, // 5 -> gprsc
    };
}

solution algorithms::run_experiment(const int algorithm, const int budget) {
    solution out;

    if (algorithm >= 0 && algorithm < static_cast<int>(algorithm_functions.size())) {
        out = algorithm_functions[algorithm](*this, budget);
    } else {
        throw invalid_argument("Invalid algorithm index: " + to_string(algorithm));
    }

    if (out.cost > budget) {
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
    if (algorithm == 4 || algorithm == 5) { // oprsc, gprsc
        string reason;
        if (!util::is_partial_row_single_column_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid partial-row single-column solution: " + reason);
        }
    }

    return out;
}

solution algorithms::opt_full_row(const int budget) const {
    solution out;
    out.algorithm_key = "ofr";

    const RewardGrid &rewards = dep.rewards();
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols();

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        return out;
    }

    // OPTFRI-style implementation under the unified cost metric:
    // choose an optimal subset of full rows; cost depends only on max selected row and number of row crossings.
    // Each full-row crossing on the full graph costs (n+1) edges (from corridor to corridor).
    const int row_cross_cost = num_internal_cols + 1;

    vector<int> row_rewards(num_rows, 0);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_internal_cols; ++j) {
            row_rewards[i] += rewards[i][j];
        }
    }

    int best_reward = 0;
    int best_cost = 0;
    vector<int> best_rows;

    for (int max_row = 0; max_row < num_rows; ++max_row) {
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
            const int candidate_cost = 2 * max_row + crossings * row_cross_cost;

            if (candidate_cost > budget) {
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

            int candidate_reward = 0;
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
    out.cycle = util::build_full_row_cycle_from_selection(out.traversed_rows, num_internal_cols);
    out.cost = util::compute_cycle_cost(out.cycle);
    return out;
}

solution algorithms::opt_partial_row_single_column(const int budget) const {
    return solve_oprsc_left_on_grid(dep.rewards(), budget, "oprsc", false);
}

solution algorithms::opt_partial_row_single_column_right(const int budget) const {
    solution out;
    out.algorithm_key = "oprsc-right";

    const RewardGrid &rewards = dep.rewards();
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols();

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        return out;
    }

    const int right_corridor = num_internal_cols + 1;

    // Suffix rewards by number of visited internal columns from the right:
    // suff[i][c] = reward collected on row i by visiting the last c internal columns, c in [0..n].
    vector<vector<int>> suff(num_rows, vector<int>(num_internal_cols + 1, 0));
    for (int i = 0; i < num_rows; ++i) {
        for (int c = 1; c <= num_internal_cols; ++c) {
            suff[i][c] = suff[i][c - 1] + rewards[i][num_internal_cols - c];
        }
    }

    const int B2 = budget / 2;
    if (B2 <= 0) {
        out.cycle = {{0, right_corridor}};
        return out;
    }

    const int NEG_INF = -1000000000;
    vector<vector<int>> R(num_rows, vector<int>(B2 + 1, NEG_INF));
    vector<vector<int>> S(num_rows, vector<int>(B2 + 1, 0)); // chosen suffix length c in [0..n]

    for (int b = 0; b <= B2; ++b) {
        const int best_c = min(num_internal_cols, b);
        R[0][b] = suff[0][best_c];
        S[0][b] = best_c;
    }

    for (int i = 1; i < num_rows; ++i) {
        for (int b = 0; b <= B2; ++b) {
            if (b < i) {
                R[i][b] = NEG_INF;
                S[i][b] = 0;
                continue;
            }

            if (b == i) {
                R[i][b] = 0;
                S[i][b] = 0;
                continue;
            }

            int best_val = NEG_INF;
            int best_c = 0;
            for (int c = 0; c <= num_internal_cols; ++c) {
                const int idx = b - c - 1;
                if (idx < 0) {
                    continue;
                }
                if (R[i - 1][idx] <= NEG_INF / 2) {
                    continue;
                }
                const int cand = R[i - 1][idx] + suff[i][c];
                if (cand > best_val) {
                    best_val = cand;
                    best_c = c;
                }
            }
            R[i][b] = best_val;
            S[i][b] = best_c;
        }
    }

    int best_max_row = -1;
    int best_reward = 0;
    for (int i = 0; i < num_rows; ++i) {
        if (R[i][B2] > best_reward) {
            best_reward = R[i][B2];
            best_max_row = i;
        }
    }

    if (best_max_row < 0 || best_reward <= 0) {
        out.cycle = {{0, right_corridor}};
        return out;
    }

    vector<int> c_by_row(best_max_row + 1, 0);
    int b = B2;
    for (int i = best_max_row; i >= 0; --i) {
        const int c = S[i][b];
        c_by_row[i] = c;
        if (i > 0) {
            b -= (c + 1);
            if (b < 0) {
                throw runtime_error("oprsc-right backtracking produced negative reduced budget");
            }
        }
    }

    vector<pair<int, int>> cycle;
    cycle.reserve(3 * (best_max_row + 1) + 2);
    cycle.push_back({0, right_corridor});
    for (int i = 0; i <= best_max_row; ++i) {
        if (cycle.back() != make_pair(i, right_corridor)) {
            cycle.push_back({i, right_corridor});
        }
        if (c_by_row[i] > 0) {
            const int endpoint_col = right_corridor - c_by_row[i];
            cycle.push_back({i, endpoint_col});
            cycle.push_back({i, right_corridor});
            out.traversed_rows.push_back(i);
        }
        out.selected_columns_per_row.push_back(c_by_row[i]);
    }
    if (cycle.back() != make_pair(0, right_corridor)) {
        cycle.push_back({0, right_corridor});
    }

    out.reward = best_reward;
    out.cycle = move(cycle);
    out.cost = util::compute_cycle_cost(out.cycle);
    return out;
}

solution algorithms::heuristic_partial_row(const int budget) const {
    return heuristic_partial_row_impl(budget, 0);
}

solution algorithms::heuristic_partial_row_s1(const int budget) const {
    return heuristic_partial_row_impl(budget, 1);
}

solution algorithms::heuristic_partial_row_s3(const int budget) const {
    return heuristic_partial_row_impl(budget, 2);
}

solution algorithms::heuristic_partial_row_impl(const int budget, const int component_mode) const {
    // Incremental paper-driven implementation:
    // - APX baseline = max(S_L, S_R, S_F)
    // - S1 (current implementation): OFR(B) + left OPRSC on residual budget over a filtered/zeroed matrix,
    //   then, if budget remains, right OPRSC on the remaining (zeroed) rewards.
    solution apx = apx_partial_row(budget);

    solution s1 = opt_full_row(budget);
    const RewardGrid &base_rewards = dep.rewards();
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols();
    int residual_budget = budget - s1.cost;

    RewardGrid residual_grid = base_rewards;
    int max_full_row = -1;
    for (const int r : s1.traversed_rows) {
        max_full_row = max(max_full_row, r);
    }

    // Remove (zero) rows already collected by OFR and rows deeper than the last OFR row.
    vector<char> is_full_row_selected(num_rows, 0);
    for (const int r : s1.traversed_rows) {
        if (r >= 0 && r < num_rows) {
            is_full_row_selected[r] = 1;
        }
    }
    for (int i = 0; i < num_rows; ++i) {
        if ((max_full_row >= 0 && i > max_full_row) || is_full_row_selected[i]) {
            fill(residual_grid[i].begin(), residual_grid[i].end(), 0.0);
        }
    }

    auto append_cycle_from_depot = [](solution &dst, const solution &src) {
        if (src.cycle.empty()) {
            return;
        }
        if (dst.cycle.empty()) {
            dst.cycle = src.cycle;
            return;
        }
        for (size_t i = 1; i < src.cycle.size(); ++i) {
            dst.cycle.push_back(src.cycle[i]);
        }
    };

    auto zero_rewards_collected_by_left_partial = [&](RewardGrid &grid_to_zero, const solution &partial_left) {
        for (size_t i = 0; i + 1 < partial_left.cycle.size(); ++i) {
            const auto [r1, c1] = partial_left.cycle[i];
            const auto [r2, c2] = partial_left.cycle[i + 1];
            if (r1 != r2) {
                continue;
            }
            if (r1 < 0 || r1 >= num_rows) {
                continue;
            }

            if (c1 == 0 && c2 > 0 && c2 <= num_internal_cols) {
                for (int j = 0; j < c2; ++j) {
                    grid_to_zero[r1][j] = 0.0;
                }
            } else if (c2 == 0 && c1 > 0 && c1 <= num_internal_cols) {
                for (int j = 0; j < c1; ++j) {
                    grid_to_zero[r1][j] = 0.0;
                }
            }
        }
    };

    auto solve_right_rooted_on_grid = [&](const RewardGrid &grid, const int B) -> solution {
        solution out;
        out.algorithm_key = "oprsc-right";
        if (B <= 0 || grid.empty() || grid[0].empty()) {
            return out;
        }

        const int right_corridor = static_cast<int>(grid[0].size()) + 1;
        const int shift_cost = 2 * right_corridor;
        if (B < shift_cost) {
            out.cycle = {{0, 0}};
            return out;
        }

        RewardGrid flipped = grid;
        for (auto &row : flipped) {
            reverse(row.begin(), row.end());
        }
        if (!flipped.empty()) {
            fill(flipped[0].begin(), flipped[0].end(), 0.0);
        }

        solution flipped_left = solve_oprsc_left_on_grid(flipped, B - shift_cost, "oprsc-left-on-flipped", true);

        vector<pair<int, int>> cycle;
        cycle.reserve(flipped_left.cycle.size() + 4);
        for (const auto &[r, c] : flipped_left.cycle) {
            cycle.push_back({r, right_corridor - c});
        }
        if (cycle.empty()) {
            cycle = {{0, 0}};
        } else {
            vector<pair<int, int>> rooted;
            rooted.push_back({0, 0});
            rooted.push_back({0, right_corridor});
            for (size_t i = 1; i < cycle.size(); ++i) {
                rooted.push_back(cycle[i]);
            }
            if (rooted.back() != make_pair(0, right_corridor)) {
                rooted.push_back({0, right_corridor});
            }
            rooted.push_back({0, 0});

            vector<pair<int, int>> cleaned;
            for (const auto &p : rooted) {
                if (cleaned.empty() || cleaned.back() != p) {
                    cleaned.push_back(p);
                }
            }
            cycle = move(cleaned);
        }

        out.reward = flipped_left.reward;
        if (!grid.empty()) {
            for (const double v : grid[0]) {
                out.reward += v; // top-row transfer collects the full first row
            }
        }
        out.cycle = move(cycle);
        out.cost = util::compute_cycle_cost(out.cycle);
        out.traversed_rows = flipped_left.traversed_rows;
        out.selected_columns_per_row = flipped_left.selected_columns_per_row;
        return out;
    };

    auto build_full_solution_from_rows = [&](vector<int> rows_0based, const string &tag) -> solution {
        solution s;
        s.algorithm_key = "hpr";
        sort(rows_0based.begin(), rows_0based.end());
        rows_0based.erase(unique(rows_0based.begin(), rows_0based.end()), rows_0based.end());
        s.traversed_rows = rows_0based;
        s.cycle = util::build_full_row_cycle_from_selection(s.traversed_rows, num_internal_cols);
        s.cost = util::compute_cycle_cost(s.cycle);
        for (const int r : s.traversed_rows) {
            if (r < 0 || r >= num_rows) {
                continue;
            }
            for (const double v : base_rewards[r]) {
                s.reward += v;
            }
        }
        s.notes.push_back(tag);
        return s;
    };

    auto build_residual_grid_after_full_rows = [&](const vector<int> &full_rows) -> RewardGrid {
        RewardGrid grid = base_rewards;
        int k = -1;
        vector<char> selected(num_rows, 0);
        for (const int r : full_rows) {
            if (r >= 0 && r < num_rows) {
                selected[r] = 1;
                k = max(k, r);
            }
        }
        for (int i = 0; i < num_rows; ++i) {
            if ((k >= 0 && i > k) || selected[i]) {
                fill(grid[i].begin(), grid[i].end(), 0.0);
            }
        }
        return grid;
    };

    auto augment_with_residual_partials = [&](solution base_full, const string &tag) -> solution {
        int b_res = budget - base_full.cost;
        if (b_res <= 0) {
            base_full.notes.push_back(tag + ": no residual budget");
            return base_full;
        }

        RewardGrid grid = build_residual_grid_after_full_rows(base_full.traversed_rows);
        solution left_residual = solve_oprsc_left_on_grid(grid, b_res, "hgc-left-residual", false);
        b_res -= left_residual.cost;
        base_full.reward += left_residual.reward;
        append_cycle_from_depot(base_full, left_residual);
        zero_rewards_collected_by_left_partial(grid, left_residual);

        if (b_res > 0) {
            solution right_residual = solve_right_rooted_on_grid(grid, b_res);
            b_res -= right_residual.cost;
            base_full.reward += right_residual.reward;
            append_cycle_from_depot(base_full, right_residual);
        }

        base_full.cost = util::compute_cycle_cost(base_full.cycle);
        base_full.notes.push_back(tag);
        base_full.notes.push_back("Residual budget after refinement=" + to_string(max(0, b_res)));
        return base_full;
    };

    auto add_left_descent_detour = [&](solution &base_sol, RewardGrid &grid, int &b_res, const string &tag) {
        if (b_res <= 1 || base_sol.cycle.size() < 2) {
            return;
        }
        if (base_sol.cycle.back() != make_pair(0, 0)) {
            return;
        }
        const auto attach = base_sol.cycle[base_sol.cycle.size() - 2];
        const int k = attach.first;
        const int attach_col = attach.second;
        if (k < 0 || k >= num_rows || attach_col != 0) {
            return; // currently implemented only for final return on left corridor
        }

        const int b2 = b_res / 2;
        if (b2 <= 0) {
            return;
        }

        vector<vector<double>> pref(k + 1, vector<double>(num_internal_cols + 1, 0.0));
        for (int i = 0; i <= k; ++i) {
            for (int c = 1; c <= num_internal_cols; ++c) {
                pref[i][c] = pref[i][c - 1] + grid[i][c - 1];
            }
        }

        // Multiple-choice knapsack on the final left-corridor descent:
        // extra cost relative to the baseline return k->0 is exactly 2*c per row prefix.
        const double NEG_INF = -1e100;
        vector<vector<double>> dp(k + 2, vector<double>(b2 + 1, NEG_INF));
        vector<vector<int>> take(k + 1, vector<int>(b2 + 1, 0));
        dp[0][0] = 0.0;

        for (int i = 0; i <= k; ++i) {
            for (int used = 0; used <= b2; ++used) {
                if (dp[i][used] <= NEG_INF / 2) {
                    continue;
                }
                for (int c = 0; c <= num_internal_cols && used + c <= b2; ++c) {
                    const double cand = dp[i][used] + pref[i][c];
                    if (cand > dp[i + 1][used + c] + 1e-12) {
                        dp[i + 1][used + c] = cand;
                        take[i][used + c] = c;
                    } else if (fabs(cand - dp[i + 1][used + c]) <= 1e-12 && c < take[i][used + c]) {
                        take[i][used + c] = c;
                    }
                }
            }
        }

        int best_used = 0;
        double best_gain = 0.0;
        for (int used = 0; used <= b2; ++used) {
            if (dp[k + 1][used] > best_gain + 1e-12) {
                best_gain = dp[k + 1][used];
                best_used = used;
            } else if (fabs(dp[k + 1][used] - best_gain) <= 1e-12 && used < best_used) {
                best_used = used;
            }
        }
        if (best_gain <= 1e-12) {
            return;
        }

        vector<int> c_by_row(k + 1, 0);
        int used = best_used;
        for (int i = k; i >= 0; --i) {
            int c = take[i][used];
            c_by_row[i] = c;
            used -= c;
        }

        // Trim useless zero tails to keep the detour tight.
        for (int i = 0; i <= k; ++i) {
            while (c_by_row[i] > 0 && fabs(grid[i][c_by_row[i] - 1]) <= 1e-12) {
                c_by_row[i]--;
            }
        }

        // Remove the final depot point, insert detours while descending on left corridor, then close at depot once.
        base_sol.cycle.pop_back();
        for (int i = k; i >= 0; --i) {
            if (base_sol.cycle.back() != make_pair(i, 0)) {
                base_sol.cycle.push_back({i, 0});
            }
            if (c_by_row[i] > 0) {
                base_sol.cycle.push_back({i, c_by_row[i]});
                base_sol.cycle.push_back({i, 0});
                for (int j = 0; j < c_by_row[i]; ++j) {
                    grid[i][j] = 0.0;
                }
            }
        }
        if (base_sol.cycle.back() != make_pair(0, 0)) {
            base_sol.cycle.push_back({0, 0});
        }

        const int extra_cost = 2 * best_used;
        b_res -= extra_cost;
        base_sol.reward += best_gain;
        base_sol.cost = util::compute_cycle_cost(base_sol.cycle);
        base_sol.notes.push_back(tag + ": left descent detour gain=" + to_string(best_gain) +
                                 ", extra_cost=" + to_string(extra_cost));
    };

    if (residual_budget > 0) {
        // First exploit the residual on the final OFR return along the left corridor (this avoids wasting slack).
        add_left_descent_detour(s1, residual_grid, residual_budget, "S1");

        solution left_residual = solve_oprsc_left_on_grid(residual_grid, residual_budget, "s1-left-residual", false);
        s1.notes.push_back("S1 left rooted residual: reward=" + to_string(left_residual.reward) +
                           ", cost=" + to_string(left_residual.cost));
        residual_budget -= left_residual.cost;
        s1.reward += left_residual.reward;
        append_cycle_from_depot(s1, left_residual);
        zero_rewards_collected_by_left_partial(residual_grid, left_residual);

        if (residual_budget > 0) {
            solution right_residual = solve_right_rooted_on_grid(residual_grid, residual_budget);
            s1.notes.push_back("S1 right rooted residual: reward=" + to_string(right_residual.reward) +
                               ", cost=" + to_string(right_residual.cost));
            residual_budget -= right_residual.cost;
            s1.reward += right_residual.reward;
            append_cycle_from_depot(s1, right_residual);
            s1.notes.push_back("S1 added right residual phase");
        }

        s1.cost = util::compute_cycle_cost(s1.cycle);
        s1.notes.push_back("S1 = OFR + OPRSC(left residual) + optional OPRSC(right residual)");
        s1.notes.push_back("Residual budget after S1=" + to_string(max(0, residual_budget)));
    } else {
        s1.notes.push_back("S1 = OFR (no residual budget)");
    }
    s1.algorithm_key = "hpr";

    // S2: best among S_L, S_R and a combined S' built from rows suggested by both sides.
    solution s_left_budget = opt_partial_row_single_column(budget);
    solution s_right_budget = opt_partial_row_single_column_right_public(budget);
    solution s2 = s_left_budget;
    s2.algorithm_key = "hpr";
    s2.notes.push_back("S2 initialized with S_L");
    if (s_right_budget.reward > s2.reward ||
        (s_right_budget.reward == s2.reward && s_right_budget.cost < s2.cost)) {
        s2 = s_right_budget;
        s2.algorithm_key = "hpr";
        s2.notes.push_back("S2 selected S_R");
    }

    // Build S' from overlap between left and right partial selections (using selected cardinalities per row).
    vector<double> row_sum(num_rows, 0.0);
    for (int i = 0; i < num_rows; ++i) {
        for (const double v : base_rewards[i]) {
            row_sum[i] += v;
        }
    }
    vector<int> full_rows_candidate;
    int threshold = max(1, num_internal_cols / 2);
    while (threshold >= 1 && full_rows_candidate.size() < 2) {
        full_rows_candidate.clear();
        for (int i = 0; i < num_rows; ++i) {
            int left_cnt = (i < static_cast<int>(s_left_budget.selected_columns_per_row.size()))
                               ? s_left_budget.selected_columns_per_row[i]
                               : 0;
            int right_cnt = (i < static_cast<int>(s_right_budget.selected_columns_per_row.size()))
                                ? s_right_budget.selected_columns_per_row[i]
                                : 0;
            if (left_cnt + right_cnt >= threshold) {
                full_rows_candidate.push_back(i);
            }
        }
        threshold--;
    }
    if (!full_rows_candidate.empty()) {
        sort(full_rows_candidate.begin(), full_rows_candidate.end(),
             [&](int a, int b) { return row_sum[a] > row_sum[b]; });
        if (full_rows_candidate.size() % 2 == 1) {
            full_rows_candidate.pop_back();
        }
        if (!full_rows_candidate.empty()) {
            solution s_prime = build_full_solution_from_rows(full_rows_candidate, "S2 prime full rows");
            if (s_prime.cost <= budget) {
                s_prime = augment_with_residual_partials(s_prime, "S2 prime + residual partials");
                s_prime.algorithm_key = "hpr";
                if (s_prime.reward > s2.reward ||
                    (s_prime.reward == s2.reward && s_prime.cost < s2.cost)) {
                    s2 = s_prime;
                    s2.notes.push_back("S2 selected S'");
                }
            }
        }
    }

    // S3: iterate prefixes, choose two best rows up to i, then refine with partial rows.
    solution s3;
    s3.algorithm_key = "hpr";
    for (int i = 0; i < num_rows; ++i) {
        vector<int> prefix_rows;
        for (int r = 0; r <= i; ++r) {
            prefix_rows.push_back(r);
        }
        sort(prefix_rows.begin(), prefix_rows.end(),
             [&](int a, int b) {
                 if (fabs(row_sum[a] - row_sum[b]) <= 1e-12) {
                     return a < b;
                 }
                 return row_sum[a] > row_sum[b];
             });
        if (prefix_rows.size() < 2) {
            continue;
        }
        vector<int> chosen = {prefix_rows[0], prefix_rows[1]};
        solution cand = build_full_solution_from_rows(chosen, "S3 base (two best rows up to i)");
        if (cand.cost > budget) {
            continue;
        }
        cand = augment_with_residual_partials(cand, "S3 + residual partials");
        cand.algorithm_key = "hpr";
        if (cand.reward > s3.reward ||
            (cand.reward == s3.reward && cand.cost < s3.cost)) {
            s3 = cand;
            s3.notes.push_back("S3 best updated at i=" + to_string(i));
        }
    }

    if (component_mode == 1) {
        s1.algorithm_key = "hpr-s1";
        return s1;
    }
    if (component_mode == 2) {
        s3.algorithm_key = "hpr-s3";
        return s3;
    }

    solution out = apx;
    out.algorithm_key = "hpr";
    out.notes.push_back("Candidate: APX");

    auto consider_candidate = [&](const solution &cand, const string &name) {
        if (cand.reward > out.reward ||
            (cand.reward == out.reward && cand.cost < out.cost)) {
            out = cand;
            out.algorithm_key = "hpr";
            out.notes.push_back("Selected candidate: " + name);
        }
    };

    consider_candidate(s1, "S1");
    consider_candidate(s3, "S3");
    out.notes.push_back("HPR currently uses APR, S1, S3 (S2 removed)");
    out.notes.push_back("TODO: refine S3 to match paper exactly (current implementation is paper-inspired)");
    return out;
}

solution algorithms::apx_partial_row(const int budget) const {
    solution out;
    out.algorithm_key = "apr";

    const int num_internal_cols = dep.cols();
    const int right_shift_cost = 2 * (num_internal_cols + 1);

    solution s_left = opt_partial_row_single_column(budget);
    solution s_full = opt_full_row(budget);

    solution s_right_rooted;
    s_right_rooted.algorithm_key = "oprsc-right-rooted";
    if (budget >= right_shift_cost) {
        s_right_rooted = opt_partial_row_single_column_right_public(budget);
        s_right_rooted.algorithm_key = "oprsc-right-rooted";
    }

    const vector<solution *> candidates = {&s_left, &s_right_rooted, &s_full};
    const vector<string> candidate_names = {"S_L(oprsc-left)", "S_R(oprsc-right)", "S_F(ofr)"};

    int best_idx = 0;
    for (int i = 1; i < static_cast<int>(candidates.size()); ++i) {
        const solution &cur = *candidates[i];
        const solution &best = *candidates[best_idx];
        if (cur.reward > best.reward) {
            best_idx = i;
        } else if (cur.reward == best.reward && cur.cost < best.cost) {
            best_idx = i;
        }
    }

    out = *candidates[best_idx];
    out.algorithm_key = "apr";
    out.notes.push_back("Selected candidate: " + candidate_names[best_idx]);
    return out;
}

solution algorithms::greedy_full_row(const int budget) const {
    solution out;
    out.algorithm_key = "gfr";

    const RewardGrid &reward_map = dep.rewards();
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols(); // internal columns (same matrix convention as legacy code)

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        return out;
    }

    vector<pair<int, int>> rows; // {row_reward, row_id_1based}
    rows.reserve(num_rows);
    for (int i = 0; i < num_rows; ++i) {
        int row_reward = 0;
        for (int j = 0; j < num_internal_cols; ++j) {
            row_reward += reward_map[i][j];
        }
        rows.emplace_back(row_reward, i + 1);
    }

    sort(rows.begin(), rows.end(),
         [](const pair<int, int> &a, const pair<int, int> &b) {
             if (a.first == b.first) {
                 return a.second < b.second;
             }
             return a.first > b.first;
         });

    int total_reward = 0;
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
            util::build_full_row_cycle_from_selection(candidate_rows_0based, num_internal_cols);
        const int candidate_cost = util::compute_cycle_cost(candidate_cycle);

        if (candidate_cost <= budget) {
            selected_rows_1based.push_back(row_id);
            total_reward += row_reward;
        }
    }

    out.reward = total_reward;
    out.traversed_rows.reserve(selected_rows_1based.size());
    for (const int row_id : selected_rows_1based) {
        out.traversed_rows.push_back(row_id - 1);
    }

    out.cycle = util::build_full_row_cycle_from_selection(out.traversed_rows, num_internal_cols);
    out.cost = util::compute_cycle_cost(out.cycle);

    return out;
}

solution algorithms::greedy_partial_row_single_column(const int budget) const {
    solution out;
    out.algorithm_key = "gprsc";

    const RewardGrid &rewards = dep.rewards();
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols(); // internal profitable columns

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        return out;
    }

    // Prefix rewards on internal columns: pref[i][c] is reward from first c internal columns (c in [0..n]).
    vector<vector<int>> pref(num_rows, vector<int>(num_internal_cols + 1, 0));
    for (int i = 0; i < num_rows; ++i) {
        for (int c = 1; c <= num_internal_cols; ++c) {
            pref[i][c] = pref[i][c - 1] + rewards[i][c - 1];
        }
    }

    int current_row = 0; // always on left corridor
    int total_cost = 0;
    int total_reward = 0;

    // selected_prefix[i] = number of internal columns already collected on row i (0..n)
    vector<int> selected_prefix(num_rows, 0);
    vector<pair<int, int>> cycle;
    cycle.push_back({0, 0});

    while (true) {
        double best_score = -1.0;
        int best_row = -1;
        int best_c = 0;

        for (int i = 0; i < num_rows; ++i) {
            for (int c = selected_prefix[i] + 1; c <= num_internal_cols; ++c) {
                const int add_reward = pref[i][c] - pref[i][selected_prefix[i]];
                if (add_reward <= 0) {
                    continue;
                }

                const int add_cost = abs(i - current_row) + 2 * c;
                const int return_home_cost = i;
                if (total_cost + add_cost + return_home_cost > budget) {
                    continue;
                }

                if (add_cost <= 0) {
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
        total_cost += abs(best_row - current_row);
        if (cycle.back() != make_pair(best_row, 0)) {
            cycle.push_back({best_row, 0});
        }

        // Partial-row excursion on the full graph: to internal column 'best_c' and back to left corridor.
        cycle.push_back({best_row, best_c});
        cycle.push_back({best_row, 0});
        total_cost += 2 * best_c;

        total_reward += pref[best_row][best_c] - pref[best_row][selected_prefix[best_row]];
        selected_prefix[best_row] = best_c;
        out.traversed_rows.push_back(best_row);
        current_row = best_row;
    }

    if (current_row != 0) {
        cycle.push_back({0, 0});
        total_cost += current_row;
    }

    out.reward = total_reward;
    out.selected_columns_per_row = selected_prefix;
    out.cycle = move(cycle);
    out.cost = util::compute_cycle_cost(out.cycle);

    // Keep reward/cost based on the explicit path/selection we built.
    // The path-cost recomputation should match total_cost; keep a note if not for debugging.
    if (out.cost != total_cost) {
        out.notes.push_back("Warning: internal cost accumulator mismatch with cycle cost");
    }

    return out;
}

solution algorithms::opt_partial_row_single_column_right_public(const int budget) const {
    solution out;
    out.algorithm_key = "oprsc-right";

    const int num_internal_cols = dep.cols();
    const int right_corridor = num_internal_cols + 1;
    const int shift_cost = 2 * right_corridor;

    if (budget < shift_cost) {
        // Not enough budget to go to the right corridor and come back.
        out.cycle = {{0, 0}};
        out.cost = 0;
        out.reward = 0;
        return out;
    }

    // Flip internal columns and solve left-version on reduced budget, disallowing partial collection on row 0.
    RewardGrid flipped = dep.rewards();
    for (auto &row : flipped) {
        reverse(row.begin(), row.end());
    }
    if (!flipped.empty()) {
        fill(flipped[0].begin(), flipped[0].end(), 0.0);
    }
    solution flipped_left = solve_oprsc_left_on_grid(flipped, budget - shift_cost, "oprsc-left-on-flipped", true);

    vector<pair<int, int>> cycle;
    cycle.reserve(flipped_left.cycle.size() + 4);
    for (const auto &[r, c] : flipped_left.cycle) {
        cycle.push_back({r, right_corridor - c});
    }
    // Add the fixed top-row transfer from depot to right corridor and back.
    if (cycle.empty()) {
        cycle = {{0, 0}};
    } else {
        vector<pair<int, int>> rooted;
        rooted.push_back({0, 0});
        if (rooted.back() != make_pair(0, right_corridor)) {
            rooted.push_back({0, right_corridor});
        }
        for (size_t i = 1; i < cycle.size(); ++i) {
            rooted.push_back(cycle[i]);
        }
        if (rooted.back() != make_pair(0, right_corridor)) {
            rooted.push_back({0, right_corridor});
        }
        rooted.push_back({0, 0});

        vector<pair<int, int>> cleaned;
        cleaned.reserve(rooted.size());
        for (const auto &p : rooted) {
            if (cleaned.empty() || cleaned.back() != p) {
                cleaned.push_back(p);
            }
        }
        cycle = move(cleaned);
    }

    out.reward = flipped_left.reward;
    if (!dep.rewards().empty()) {
        int top_row_reward = 0;
        for (const double value : dep.rewards()[0]) {
            top_row_reward += value;
        }
        out.reward += top_row_reward;
    }
    out.cycle = move(cycle);
    out.cost = util::compute_cycle_cost(out.cycle);
    out.traversed_rows = flipped_left.traversed_rows;
    out.selected_columns_per_row = flipped_left.selected_columns_per_row;
    return out;
}
