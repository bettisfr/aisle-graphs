#include "algorithms.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../util/util.h"

#include <gurobi_c++.h>

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
        &algorithms::opt_partial_row,               // 1 -> opr
        &algorithms::greedy_full_row,               // 2 -> gfr

        &algorithms::heuristic_partial_row,         // 3 -> hpr
        &algorithms::apx_partial_row,               // 4 -> apr

        &algorithms::opt_partial_row_single_column, // 5 -> oprsc
        &algorithms::greedy_partial_row_single_column, // 6 -> gprsc
        &algorithms::greedy_partial_row,            // 7 -> gpr
        &algorithms::opt_partial_row_dp,            // 8 -> opr-dp
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
    if (algorithm == 0 || algorithm == 2) { // ofr, gfr
        string reason;
        if (!util::is_full_row_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid full-row solution: " + reason);
        }
    }
    if (algorithm == 5 || algorithm == 6) { // oprsc, gprsc
        string reason;
        if (!util::is_partial_row_single_column_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid partial-row single-column solution: " + reason);
        }
    }
    if (algorithm == 1 || algorithm == 3 || algorithm == 4 || algorithm == 7 || algorithm == 8) { // opr, hpr, apr, gpr, opr-dp
        string reason;
        if (!util::is_partial_row_solution_feasible(out, dep.rows(), dep.cols(), budget, &reason)) {
            throw runtime_error("Invalid aisle-graph path solution: " + reason);
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

solution algorithms::solve_opr_ilp_gurobi(const int budget) const {
    solution out;
    out.algorithm_key = "opr";
    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols();
    const int num_cols = num_internal_cols + 2; // include left/right corridors

    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        out.cycle = {{0, 0}};
        return out;
    }

    vector<vector<int>> full_reward(num_rows, vector<int>(num_cols, 0));
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_internal_cols; ++j) {
            full_reward[i][j + 1] = static_cast<int>(dep.rewards()[i][j]);
        }
    }

    struct Arc {
        int from;
        int to;
    };

    const int num_nodes = num_rows * num_cols;
    const int depot = 0; // (0,0)
    auto node_id = [num_cols](const int r, const int c) { return r * num_cols + c; };
    auto node_rc = [num_cols](const int v) { return pair<int, int>{v / num_cols, v % num_cols}; };

    // Aisle-graph shortest distance from depot (0,0) to node v.
    // Vertical moves are only on external corridors.
    vector<int> dist_from_depot(num_nodes, 0);
    for (int v = 0; v < num_nodes; ++v) {
        const auto [r, c] = node_rc(v);
        const int via_left = r + c;
        const int via_right = (num_cols - 1) + r + (num_cols - 1 - c);
        dist_from_depot[v] = min(via_left, via_right);
    }

    vector<Arc> arcs;
    arcs.reserve(4 * num_nodes);
    vector<vector<int>> out_arcs(num_nodes), in_arcs(num_nodes);
    for (int r = 0; r < num_rows; ++r) {
        for (int c = 0; c < num_cols; ++c) {
            const int u = node_id(r, c);

            // Horizontal edges are allowed on every row.
            if (c - 1 >= 0) {
                const int v = node_id(r, c - 1);
                const int aidx = static_cast<int>(arcs.size());
                arcs.push_back({u, v});
                out_arcs[u].push_back(aidx);
                in_arcs[v].push_back(aidx);
            }
            if (c + 1 < num_cols) {
                const int v = node_id(r, c + 1);
                const int aidx = static_cast<int>(arcs.size());
                arcs.push_back({u, v});
                out_arcs[u].push_back(aidx);
                in_arcs[v].push_back(aidx);
            }

            // Vertical edges are allowed only on external corridors.
            const bool external_col = (c == 0 || c == num_cols - 1);
            if (external_col) {
                if (r - 1 >= 0) {
                    const int v = node_id(r - 1, c);
                    const int aidx = static_cast<int>(arcs.size());
                    arcs.push_back({u, v});
                    out_arcs[u].push_back(aidx);
                    in_arcs[v].push_back(aidx);
                }
                if (r + 1 < num_rows) {
                    const int v = node_id(r + 1, c);
                    const int aidx = static_cast<int>(arcs.size());
                    arcs.push_back({u, v});
                    out_arcs[u].push_back(aidx);
                    in_arcs[v].push_back(aidx);
                }
            }
        }
    }

    const int horizon = budget;
    vector<vector<char>> feasible_state(horizon + 1, vector<char>(num_nodes, 0));
    for (int t = 0; t <= horizon; ++t) {
        for (int v = 0; v < num_nodes; ++v) {
            // State (t,v) is feasible iff v is reachable by time t and can still return by horizon.
            if (dist_from_depot[v] <= t && dist_from_depot[v] <= horizon - t) {
                feasible_state[t][v] = 1;
            }
        }
    }

    try {
        GRBEnv env = GRBEnv(true);
        env.set("OutputFlag", "0");
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 300.0);

        vector<vector<GRBVar>> pos(horizon + 1, vector<GRBVar>(num_nodes));
        for (int t = 0; t <= horizon; ++t) {
            for (int v = 0; v < num_nodes; ++v) {
                const double ub = feasible_state[t][v] ? 1.0 : 0.0;
                pos[t][v] = model.addVar(0.0, ub, 0.0, GRB_BINARY, "p_" + to_string(t) + "_" + to_string(v));
            }
        }

        vector<vector<GRBVar>> move(horizon, vector<GRBVar>(arcs.size()));
        vector<GRBVar> idle(horizon);
        for (int t = 0; t < horizon; ++t) {
            idle[t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "idle_" + to_string(t));
            for (size_t a = 0; a < arcs.size(); ++a) {
                const int u = arcs[a].from;
                const int v = arcs[a].to;
                const double ub = (feasible_state[t][u] && feasible_state[t + 1][v]) ? 1.0 : 0.0;
                move[t][a] = model.addVar(0.0, ub, 0.0, GRB_BINARY,
                                          "x_" + to_string(t) + "_" + to_string(a));
            }
        }

        vector<GRBVar> visit(num_nodes);
        for (int v = 0; v < num_nodes; ++v) {
            visit[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + to_string(v));
        }

        // Initial position at depot.
        for (int v = 0; v < num_nodes; ++v) {
            if (v == depot) {
                model.addConstr(pos[0][v] == 1.0, "start_depot");
            } else if (feasible_state[0][v]) {
                model.addConstr(pos[0][v] == 0.0, "start_zero_" + to_string(v));
            }
        }

        // Exactly one position at each time step.
        for (int t = 0; t <= horizon; ++t) {
            GRBLinExpr sum_pos = 0.0;
            for (int v = 0; v < num_nodes; ++v) {
                if (feasible_state[t][v]) {
                    sum_pos += pos[t][v];
                }
            }
            model.addConstr(sum_pos == 1.0, "one_pos_" + to_string(t));
        }

        for (int t = 0; t < horizon; ++t) {
            // Outgoing action from current node.
            for (int v = 0; v < num_nodes; ++v) {
                if (!feasible_state[t][v]) {
                    continue;
                }
                GRBLinExpr outflow = 0.0;
                for (const int aidx : out_arcs[v]) {
                    if (!feasible_state[t + 1][arcs[aidx].to]) {
                        continue;
                    }
                    outflow += move[t][aidx];
                }
                if (v == depot) {
                    outflow += idle[t];
                }
                model.addConstr(outflow == pos[t][v], "flow_out_" + to_string(t) + "_" + to_string(v));
            }

            // Incoming action defines next node.
            for (int v = 0; v < num_nodes; ++v) {
                if (!feasible_state[t + 1][v]) {
                    continue;
                }
                GRBLinExpr inflow = 0.0;
                for (const int aidx : in_arcs[v]) {
                    if (!feasible_state[t][arcs[aidx].from]) {
                        continue;
                    }
                    inflow += move[t][aidx];
                }
                if (v == depot) {
                    inflow += idle[t];
                }
                model.addConstr(inflow == pos[t + 1][v], "flow_in_" + to_string(t) + "_" + to_string(v));
            }
        }

        // Must end at depot after at most budget moves (idle is allowed only at depot).
        model.addConstr(pos[horizon][depot] == 1.0, "end_depot");

        // Visit linking.
        for (int v = 0; v < num_nodes; ++v) {
            GRBLinExpr seen = 0.0;
            for (int t = 0; t <= horizon; ++t) {
                if (feasible_state[t][v]) {
                    seen += pos[t][v];
                }
            }
            model.addConstr(visit[v] <= seen, "visit_link_" + to_string(v));
        }

        GRBLinExpr obj = 0.0;
        for (int r = 0; r < num_rows; ++r) {
            for (int c = 0; c < num_cols; ++c) {
                const int v = node_id(r, c);
                obj += static_cast<double>(full_reward[r][c]) * visit[v];
            }
        }
        model.setObjective(obj);
        model.optimize();

        const int status = model.get(GRB_IntAttr_Status);
        if (status != GRB_OPTIMAL) {
            throw runtime_error("Gurobi did not return an optimal solution for OPR ILP");
        }

        vector<int> path_nodes;
        path_nodes.reserve(horizon + 1);
        for (int t = 0; t <= horizon; ++t) {
            int at = -1;
            for (int v = 0; v < num_nodes; ++v) {
                if (pos[t][v].get(GRB_DoubleAttr_X) > 0.5) {
                    at = v;
                    break;
                }
            }
            if (at < 0) {
                throw runtime_error("OPR ILP path reconstruction failed: no active node at a time step");
            }
            path_nodes.push_back(at);
        }

        vector<pair<int, int>> cycle;
        cycle.reserve(path_nodes.size());
        for (const int v : path_nodes) {
            cycle.push_back({v / num_cols, v % num_cols});
        }

        vector<pair<int, int>> cleaned;
        cleaned.reserve(cycle.size());
        for (const auto &p : cycle) {
            if (cleaned.empty() || cleaned.back() != p) {
                cleaned.push_back(p);
            }
        }
        if (cleaned.empty()) {
            cleaned.push_back({0, 0});
        }
        if (cleaned.front() != make_pair(0, 0)) {
            cleaned.insert(cleaned.begin(), {0, 0});
        }
        if (cleaned.back() != make_pair(0, 0)) {
            cleaned.push_back({0, 0});
        }

        out.cycle = std::move(cleaned);
        out.cost = util::compute_cycle_cost(out.cycle);
        out.reward = static_cast<int>(llround(model.get(GRB_DoubleAttr_ObjVal)));
        out.notes.push_back("OPR ILP via Gurobi (time-expanded exact model)");
    } catch (const GRBException &e) {
        throw runtime_error("Gurobi error in OPR ILP: code=" + to_string(e.getErrorCode()) + " msg=" + e.getMessage());
    }
    return out;
}

solution algorithms::opt_partial_row(const int budget) const {
    return solve_opr_ilp_gurobi(budget);
}

solution algorithms::opt_partial_row_dp(const int budget) const {
    solution out;
    out.algorithm_key = "opr-dp";
    out.cycle = {{0, 0}};
    out.notes.push_back("TODO: implement opt_partial_row_dp");
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

solution algorithms::heuristic_partial_row_s2(const int budget) const {
    return heuristic_partial_row_impl(budget, 2);
}

solution algorithms::heuristic_partial_row_s3(const int budget) const {
    return heuristic_partial_row_impl(budget, 3);
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

    auto better = [](const solution &a, const solution &b) {
        if (a.reward != b.reward) {
            return a.reward > b.reward;
        }
        return a.cost < b.cost;
    };

    // S1 improvement (still paper-like): try several OFR anchors and always refine with residual partials.
    // This reduces sensitivity to abrupt OFR changes between consecutive budgets.
    solution best_s1;
    best_s1.algorithm_key = "hpr";
    bool has_s1 = false;
    const int max_delta = min(80, budget);
    for (int delta = 0; delta <= max_delta; delta += 2) {
        const int anchor_budget = budget - delta;
        if (anchor_budget < 0) {
            continue;
        }
        solution base = opt_full_row(anchor_budget);
        base.algorithm_key = "hpr";
        solution cand = augment_with_residual_partials(
            base, "S1 = OFR(B-" + to_string(delta) + ") + OPRSC residuals"
        );
        cand.algorithm_key = "hpr";
        if (!has_s1 || better(cand, best_s1)) {
            best_s1 = cand;
            has_s1 = true;
        }
    }
    s1 = best_s1;

    // S2 (paper-aligned): best among S_L, S_R, and S' built from combined left/right partial picks.
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

    // Build S' from rows where combined left/right picks meet a threshold.
    // Start from ceil(n/2), then fallback to ceil(n/3), ceil(n/4), ... until at least two rows exist.
    vector<double> row_sum(num_rows, 0.0);
    for (int i = 0; i < num_rows; ++i) {
        for (const double v : base_rewards[i]) {
            row_sum[i] += v;
        }
    }

    auto collect_rows_for_threshold = [&](const int threshold) {
        vector<int> rows;
        for (int i = 0; i < num_rows; ++i) {
            int left_cnt = (i < static_cast<int>(s_left_budget.selected_columns_per_row.size()))
                               ? s_left_budget.selected_columns_per_row[i]
                               : 0;
            int right_cnt = (i < static_cast<int>(s_right_budget.selected_columns_per_row.size()))
                                ? s_right_budget.selected_columns_per_row[i]
                                : 0;
            if (left_cnt + right_cnt >= threshold) {
                rows.push_back(i);
            }
        }
        return rows;
    };

    // Multi-candidate S2': try several thresholds and row-order criteria, keep best feasible S'.
    auto try_update_s2_prime = [&](const vector<int> &rows_in, const int threshold, const string &variant_tag) {
        if (rows_in.size() < 2) {
            return;
        }

        const int max_rows_considered = min(static_cast<int>(rows_in.size()), 16);
        for (int take = 2; take <= max_rows_considered; take += 2) {
            vector<int> rows(rows_in.begin(), rows_in.begin() + take);

            solution s_prime = build_full_solution_from_rows(
                rows,
                "S2 prime full rows (threshold=" + to_string(threshold) +
                ", variant=" + variant_tag +
                ", take=" + to_string(take) + ")"
            );
            if (s_prime.cost > budget) {
                continue;
            }
            s_prime = augment_with_residual_partials(s_prime, "S2 prime + residual partials");
            s_prime.algorithm_key = "hpr";
            if (better(s_prime, s2)) {
                s2 = s_prime;
                s2.notes.push_back("S2 selected S' (" + variant_tag + ", take=" + to_string(take) + ")");
            }
        }
    };

    for (int d = 2; d <= num_internal_cols; ++d) {
        const int threshold = (num_internal_cols + d - 1) / d; // ceil(n/d)
        vector<int> cand = collect_rows_for_threshold(threshold);
        if (cand.size() < 2) {
            continue;
        }

        vector<int> by_reward = cand;
        sort(by_reward.begin(), by_reward.end(),
             [&](const int a, const int b) {
                 if (fabs(row_sum[a] - row_sum[b]) <= 1e-12) {
                     return a < b;
                 }
                 return row_sum[a] > row_sum[b];
             });
        try_update_s2_prime(by_reward, threshold, "reward");

        vector<int> by_depth = cand;
        sort(by_depth.begin(), by_depth.end(),
             [](const int a, const int b) { return a < b; });
        try_update_s2_prime(by_depth, threshold, "depth");

        vector<int> by_density = cand;
        sort(by_density.begin(), by_density.end(),
             [&](const int a, const int b) {
                 int left_a = (a < static_cast<int>(s_left_budget.selected_columns_per_row.size()))
                                  ? s_left_budget.selected_columns_per_row[a]
                                  : 0;
                 int right_a = (a < static_cast<int>(s_right_budget.selected_columns_per_row.size()))
                                   ? s_right_budget.selected_columns_per_row[a]
                                   : 0;
                 int left_b = (b < static_cast<int>(s_left_budget.selected_columns_per_row.size()))
                                  ? s_left_budget.selected_columns_per_row[b]
                                  : 0;
                 int right_b = (b < static_cast<int>(s_right_budget.selected_columns_per_row.size()))
                                   ? s_right_budget.selected_columns_per_row[b]
                                   : 0;
                 const int sum_a = left_a + right_a;
                 const int sum_b = left_b + right_b;
                 if (sum_a == sum_b) {
                     if (fabs(row_sum[a] - row_sum[b]) <= 1e-12) {
                         return a < b;
                     }
                     return row_sum[a] > row_sum[b];
                 }
                 return sum_a > sum_b;
             });
        try_update_s2_prime(by_density, threshold, "density");
    }

    // S3: iterate prefixes, evaluate multiple promising row pairs, then refine with partial rows.
    solution s3;
    s3.algorithm_key = "hpr";
    const int s3_pair_pool_size = 8;
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

        const int pool_size = min(static_cast<int>(prefix_rows.size()), s3_pair_pool_size);
        solution best_for_i;
        best_for_i.algorithm_key = "hpr";
        bool found_for_i = false;
        for (int a = 0; a < pool_size; ++a) {
            for (int b = a + 1; b < pool_size; ++b) {
                vector<int> chosen = {prefix_rows[a], prefix_rows[b]};
                solution cand = build_full_solution_from_rows(chosen, "S3 base (row-pair in top pool up to i)");
                if (cand.cost > budget) {
                    continue;
                }
                cand = augment_with_residual_partials(cand, "S3 + residual partials");
                cand.algorithm_key = "hpr";
                if (!found_for_i ||
                    cand.reward > best_for_i.reward ||
                    (cand.reward == best_for_i.reward && cand.cost < best_for_i.cost)) {
                    best_for_i = cand;
                    found_for_i = true;
                }
            }
        }
        if (!found_for_i) {
            continue;
        }
        if (best_for_i.reward > s3.reward ||
            (best_for_i.reward == s3.reward && best_for_i.cost < s3.cost)) {
            s3 = best_for_i;
            s3.notes.push_back("S3 best updated at i=" + to_string(i));
        }
    }

    if (component_mode == 1) {
        s1.algorithm_key = "hpr-s1";
        return s1;
    }
    if (component_mode == 2) {
        s2.algorithm_key = "hpr-s2";
        return s2;
    }
    if (component_mode == 3) {
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
    consider_candidate(s2, "S2");
    consider_candidate(s3, "S3");
    out.notes.push_back("HPR currently uses APR, S1, S2, S3");
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

solution algorithms::greedy_partial_row(const int budget) const {
    solution out;
    out.algorithm_key = "gpr";

    const int num_rows = dep.rows();
    const int num_internal_cols = dep.cols();
    if (budget <= 0 || num_rows <= 0 || num_internal_cols <= 0) {
        out.cycle = {{0, 0}};
        return out;
    }

    const int num_cols = num_internal_cols + 2; // includes left/right corridors
    const int ending_row = 0;

    vector<vector<int>> reward_map(num_rows, vector<int>(num_cols, 0));
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_internal_cols; ++j) {
            reward_map[i][j + 1] = static_cast<int>(dep.rewards()[i][j]);
        }
    }
    const vector<vector<int>> original_rewards = reward_map;

    int current_row = 0;
    int current_side = 1; // 1 = left, 2 = right
    int total_cost = 0;
    int total_reward = reward_map[current_row][0];
    reward_map[current_row][0] = 0;

    vector<int> tour;
    tour.push_back(1); // MATLAB-style 1-based linearized index

    vector<vector<int>> cumulative_reward_1(num_rows, vector<int>(num_cols, 0));
    vector<vector<int>> cumulative_reward_2(num_rows, vector<int>(num_cols, 0));
    vector<vector<int>> cumulative_cost_1(num_rows, vector<int>(num_cols, 0));
    vector<vector<int>> cumulative_cost_2(num_rows, vector<int>(num_cols, 0));

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            int left_sum = 0;
            for (int t = 0; t <= j; ++t) {
                left_sum += reward_map[i][t];
            }
            cumulative_reward_1[i][j] = left_sum;

            int rev_j = num_cols - 1 - j;
            int right_sum = 0;
            for (int t = num_cols - 1; t >= rev_j; --t) {
                right_sum += reward_map[i][t];
            }
            cumulative_reward_2[i][rev_j] = right_sum;

            cumulative_cost_1[i][j] = j;
            cumulative_cost_2[i][rev_j] = j;
        }
    }

    auto sum_all = [](const vector<vector<int>> &mat) {
        long long s = 0;
        for (const auto &row : mat) {
            for (const int v : row) {
                s += v;
            }
        }
        return s;
    };

    auto push_linear = [&](const int row, const int col) {
        tour.push_back(row * num_cols + col + 1);
    };

    while (sum_all(cumulative_reward_1) > 0 || sum_all(cumulative_reward_2) > 0) {
        const int budget_left = budget - total_cost;
        if (budget_left <= 0) {
            break;
        }

        if (current_side == 1) {
            double best_loop_heuristic = -1.0;
            int best_loop_index = 0;
            double best_muted_heuristic = -1.0;
            int best_muted_row = 0;
            int best_muted_col = 0;

            for (int i = 0; i < num_rows; ++i) {
                const int distance_side = abs(current_row - i);

                const int loop_den = cumulative_cost_1[i][num_cols - 1] + distance_side;
                if (loop_den > 0) {
                    const double h = static_cast<double>(cumulative_reward_1[i][num_cols - 1]) /
                                     static_cast<double>(loop_den);
                    if (h > best_loop_heuristic) {
                        best_loop_heuristic = h;
                        best_loop_index = i;
                    }
                }

                for (int j = 0; j < num_cols; ++j) {
                    const int den = 2 * cumulative_cost_1[i][j] + distance_side;
                    if (den <= 0) {
                        continue;
                    }
                    const double h = static_cast<double>(cumulative_reward_1[i][j]) / static_cast<double>(den);
                    if (h > best_muted_heuristic) {
                        best_muted_heuristic = h;
                        best_muted_row = i;
                        best_muted_col = j;
                    }
                }
            }

            const int loop_cost = abs(current_row - best_loop_index) +
                                  2 * cumulative_cost_1[best_loop_index][num_cols - 1] +
                                  abs(best_loop_index - ending_row);
            const int muted_loop_cost = abs(current_row - best_muted_row) +
                                        2 * cumulative_cost_1[best_muted_row][best_muted_col] +
                                        abs(best_muted_row - ending_row);

            if ((best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= budget_left)) {
                total_cost += abs(current_row - best_loop_index) + cumulative_cost_1[best_loop_index][num_cols - 1];
                total_reward += cumulative_reward_1[best_loop_index][num_cols - 1];

                const int step = (current_row <= best_loop_index) ? 1 : -1;
                for (int i = current_row; i != best_loop_index; i += step) {
                    push_linear(i, 0);
                    total_reward += cumulative_reward_1[i][0];
                    for (int j = 0; j < num_cols; ++j) {
                        cumulative_reward_1[i][j] = max(0, cumulative_reward_1[i][j] - cumulative_reward_1[i][0]);
                    }
                    if (num_cols > 1) {
                        cumulative_reward_2[i][0] = cumulative_reward_2[i][1];
                    }
                }

                for (int j = 0; j < num_cols; ++j) {
                    push_linear(best_loop_index, j);
                }

                current_row = best_loop_index;
                current_side = 2;
                fill(cumulative_reward_1[current_row].begin(), cumulative_reward_1[current_row].end(), 0);
                fill(cumulative_reward_2[current_row].begin(), cumulative_reward_2[current_row].end(), 0);
            } else if (muted_loop_cost <= budget_left) {
                total_cost += abs(current_row - best_muted_row) + 2 * cumulative_cost_1[best_muted_row][best_muted_col];

                const int step = (current_row <= best_muted_row) ? 1 : -1;
                for (int i = current_row; i != best_muted_row; i += step) {
                    push_linear(i, 0);
                    total_reward += cumulative_reward_1[i][0];
                    for (int j = 0; j < num_cols; ++j) {
                        cumulative_reward_1[i][j] = max(0, cumulative_reward_1[i][j] - cumulative_reward_1[i][0]);
                    }
                    if (num_cols > 1) {
                        cumulative_reward_2[i][0] = cumulative_reward_2[i][1];
                    }
                }

                for (int j = 0; j <= best_muted_col; ++j) {
                    push_linear(best_muted_row, j);
                }
                for (int j = best_muted_col; j >= 0; --j) {
                    push_linear(best_muted_row, j);
                }

                const int accumulated_reward = cumulative_reward_1[best_muted_row][best_muted_col];
                total_reward += accumulated_reward;
                current_row = best_muted_row;

                for (int j = 0; j < num_cols; ++j) {
                    cumulative_reward_1[best_muted_row][j] = max(0, cumulative_reward_1[best_muted_row][j] - accumulated_reward);
                }
                const int suffix_val = (best_muted_col + 1 < num_cols) ? cumulative_reward_2[best_muted_row][best_muted_col + 1] : 0;
                for (int j = 0; j <= best_muted_col; ++j) {
                    cumulative_reward_2[best_muted_row][j] = suffix_val;
                }
            } else {
                vector<vector<int>> reachable(num_rows, vector<int>(num_cols, 0));
                for (int i = 0; i < num_rows; ++i) {
                    const int tmp_sum = abs(current_row - i) + abs(i - ending_row);
                    for (int j = 0; j < num_cols; ++j) {
                        const int temp_cost = tmp_sum + 2 * cumulative_cost_1[i][j];
                        reachable[i][j] = (temp_cost <= budget_left) ? 1 : 0;
                    }
                }

                for (int i = 0; i < num_rows; ++i) {
                    bool any_unreachable = false;
                    for (int j = 0; j < num_cols; ++j) {
                        if (!reachable[i][j]) {
                            cumulative_reward_1[i][j] = 0;
                            any_unreachable = true;
                        }
                    }
                    if (any_unreachable) {
                        fill(cumulative_reward_2[i].begin(), cumulative_reward_2[i].end(), 0);
                    }
                }
            }
        } else { // current_side == 2
            double best_loop_heuristic = -1.0;
            int best_loop_index = 0;
            double best_muted_heuristic = -1.0;
            int best_muted_row = 0;
            int best_muted_col = 0;

            for (int i = 0; i < num_rows; ++i) {
                const int distance_side = abs(current_row - i);

                const int loop_den = cumulative_cost_2[i][0] + distance_side;
                if (loop_den > 0) {
                    const double h = static_cast<double>(cumulative_reward_2[i][0]) / static_cast<double>(loop_den);
                    if (h > best_loop_heuristic) {
                        best_loop_heuristic = h;
                        best_loop_index = i;
                    }
                }

                for (int j = 0; j < num_cols; ++j) {
                    const int den = 2 * cumulative_cost_2[i][j] + distance_side;
                    if (den <= 0) {
                        continue;
                    }
                    const double h = static_cast<double>(cumulative_reward_2[i][j]) / static_cast<double>(den);
                    if (h > best_muted_heuristic) {
                        best_muted_heuristic = h;
                        best_muted_row = i;
                        best_muted_col = j;
                    }
                }
            }

            const int loop_cost = abs(current_row - best_loop_index) +
                                  cumulative_cost_2[best_loop_index][0] +
                                  abs(best_loop_index - ending_row);
            const int muted_loop_cost = abs(current_row - best_muted_row) +
                                        2 * cumulative_cost_2[best_muted_row][best_muted_col] +
                                        abs(best_muted_row - ending_row) +
                                        cumulative_cost_2[ending_row][0];

            if ((best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= budget_left)) {
                total_cost += abs(current_row - best_loop_index) + cumulative_cost_2[best_loop_index][0];
                total_reward += cumulative_reward_2[best_loop_index][0];

                const int step = (current_row <= best_loop_index) ? 1 : -1;
                for (int i = current_row; i != best_loop_index; i += step) {
                    push_linear(i, num_cols - 1);
                    total_reward += cumulative_reward_2[i][num_cols - 1];
                    for (int j = 0; j < num_cols; ++j) {
                        cumulative_reward_2[i][j] = max(0, cumulative_reward_2[i][j] - cumulative_reward_2[i][num_cols - 1]);
                    }
                    if (num_cols > 1) {
                        cumulative_reward_1[i][num_cols - 1] = cumulative_reward_1[i][num_cols - 2];
                    }
                }

                for (int j = num_cols - 1; j >= 0; --j) {
                    push_linear(best_loop_index, j);
                }

                current_row = best_loop_index;
                current_side = 1;
                fill(cumulative_reward_1[current_row].begin(), cumulative_reward_1[current_row].end(), 0);
                fill(cumulative_reward_2[current_row].begin(), cumulative_reward_2[current_row].end(), 0);
            } else if (muted_loop_cost <= budget_left) {
                total_cost += abs(current_row - best_muted_row) + 2 * cumulative_cost_2[best_muted_row][best_muted_col];

                const int step = (current_row <= best_muted_row) ? 1 : -1;
                for (int i = current_row; i != best_muted_row; i += step) {
                    push_linear(i, num_cols - 1);
                    total_reward += cumulative_reward_2[i][num_cols - 1];
                    for (int j = 0; j < num_cols; ++j) {
                        cumulative_reward_2[i][j] = max(0, cumulative_reward_2[i][j] - cumulative_reward_2[i][num_cols - 1]);
                    }
                    if (num_cols > 1) {
                        cumulative_reward_1[i][num_cols - 1] = cumulative_reward_1[i][num_cols - 2];
                    }
                }

                for (int j = num_cols - 1; j >= best_muted_col; --j) {
                    push_linear(best_muted_row, j);
                }
                for (int j = best_muted_col; j < num_cols; ++j) {
                    push_linear(best_muted_row, j);
                }

                const int accumulated_reward = cumulative_reward_2[best_muted_row][best_muted_col];
                total_reward += accumulated_reward;
                current_row = best_muted_row;

                for (int j = 0; j < num_cols; ++j) {
                    cumulative_reward_2[best_muted_row][j] = max(0, cumulative_reward_2[best_muted_row][j] - accumulated_reward);
                }
                const int prefix_val = (best_muted_col - 1 >= 0) ? cumulative_reward_1[best_muted_row][best_muted_col - 1] : 0;
                for (int j = num_cols - 1; j >= best_muted_col; --j) {
                    cumulative_reward_1[best_muted_row][j] = prefix_val;
                }
            } else {
                vector<int> reachable(num_rows, 0);
                for (int i = 0; i < num_rows; ++i) {
                    const int temp_cost = abs(current_row - i) + cumulative_cost_2[i][0] + abs(i - ending_row);
                    reachable[i] = (temp_cost <= budget_left) ? 1 : 0;
                }
                for (int i = 0; i < num_rows; ++i) {
                    if (!reachable[i]) {
                        fill(cumulative_reward_1[i].begin(), cumulative_reward_1[i].end(), 0);
                        fill(cumulative_reward_2[i].begin(), cumulative_reward_2[i].end(), 0);
                    }
                }
                for (int j = 0; j <= best_muted_col; ++j) {
                    cumulative_reward_2[best_muted_row][j] = 0;
                }
                const int drop = cumulative_reward_1[best_muted_row][best_muted_col];
                for (int j = 0; j < num_cols; ++j) {
                    cumulative_reward_1[best_muted_row][j] = max(0, cumulative_reward_1[best_muted_row][j] - drop);
                }
            }
        }
    }

    if (current_side == 1) {
        total_cost += abs(current_row - ending_row);
        const int step = (current_row > ending_row) ? -1 : 1;
        for (int i = current_row;; i += step) {
            push_linear(i, 0);
            if (i == ending_row) {
                break;
            }
        }
    } else {
        total_cost += abs(current_row - ending_row) + cumulative_cost_2[current_row][0];
        for (int j = num_cols - 1; j >= 0; --j) {
            push_linear(current_row, j);
        }
        const int step = (current_row > ending_row) ? -1 : 1;
        for (int i = current_row;; i += step) {
            push_linear(i, 0);
            if (i == ending_row) {
                break;
            }
        }
    }

    vector<int> dedup_tour;
    dedup_tour.reserve(tour.size());
    for (const int idx : tour) {
        if (dedup_tour.empty() || dedup_tour.back() != idx) {
            dedup_tour.push_back(idx);
        }
    }

    out.cycle.clear();
    out.cycle.reserve(dedup_tour.size());
    for (const int idx1 : dedup_tour) {
        const int zero = idx1 - 1;
        const int row = zero / num_cols;
        const int col = zero % num_cols;
        out.cycle.push_back({row, col});
    }
    if (out.cycle.empty()) {
        out.cycle.push_back({0, 0});
    }

    vector<vector<char>> picked(num_rows, vector<char>(num_cols, 0));
    int reward_from_cycle = 0;
    for (const auto &[r, c] : out.cycle) {
        if (r < 0 || r >= num_rows || c < 0 || c >= num_cols) {
            continue;
        }
        if (!picked[r][c]) {
            picked[r][c] = 1;
            reward_from_cycle += original_rewards[r][c];
        }
    }

    out.reward = reward_from_cycle;
    out.cost = util::compute_cycle_cost(out.cycle);
    if (out.cost != total_cost) {
        out.notes.push_back("Warning: gpr cost accumulator mismatch with reconstructed cycle cost");
    }
    if (out.reward != total_reward) {
        out.notes.push_back("Info: gpr reward recomputed from cycle differs from accumulator");
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
