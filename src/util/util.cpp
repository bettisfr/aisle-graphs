#include <algorithm>
#include <cmath>
#include <set>

#include "util.h"

int util::compute_cycle_cost(const vector<pair<int, int>> &cycle) {
    if (cycle.empty()) {
        return 0;
    }

    long long total = 0;
    for (size_t i = 0; i + 1 < cycle.size(); ++i) {
        const auto [r1, c1] = cycle[i];
        const auto [r2, c2] = cycle[i + 1];
        total += llabs(static_cast<long long>(r2) - static_cast<long long>(r1));
        total += llabs(static_cast<long long>(c2) - static_cast<long long>(c1));
    }
    return static_cast<int>(total);
}

vector<pair<int, int>> util::build_full_row_cycle_from_selection(const vector<int> &selected_rows_0based,
                                                                 const int internal_cols) {
    vector<pair<int, int>> cycle;
    if (internal_cols <= 0) {
        return cycle;
    }

    const int left_corridor = 0;
    const int right_corridor = internal_cols + 1;
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

bool util::is_full_row_solution_feasible(const solution &sol, const int rows, const int internal_cols, const int budget,
                                         string *reason) {
    auto fail = [&](const string &msg) {
        if (reason != nullptr) {
            *reason = msg;
        }
        return false;
    };

    if (rows <= 0 || internal_cols <= 0) {
        return fail("Invalid instance dimensions");
    }

    const int left_corridor = 0;
    const int right_corridor = internal_cols + 1;
    const int graph_cols = internal_cols + 2;

    if (sol.cost > budget) {
        return fail("Solution cost exceeds budget");
    }

    if (sol.cycle.empty()) {
        // Empty solution is feasible only if it also reports zero cost and no traversed rows.
        if (sol.cost != 0 || !sol.traversed_rows.empty()) {
            return fail("Empty cycle with non-empty rows or non-zero cost");
        }
        return true;
    }

    if (sol.cycle.front() != make_pair(0, left_corridor) || sol.cycle.back() != make_pair(0, left_corridor)) {
        return fail("Cycle must start/end at depot (0,0)");
    }

    set<int> rows_crossed_horizontally;

    for (size_t i = 0; i < sol.cycle.size(); ++i) {
        const auto [r, c] = sol.cycle[i];
        if (r < 0 || r >= rows) {
            return fail("Cycle row out of bounds");
        }
        if (c < 0 || c >= graph_cols) {
            return fail("Cycle column out of bounds");
        }

        if (i + 1 == sol.cycle.size()) {
            continue;
        }

        const auto [r2, c2] = sol.cycle[i + 1];
        const int dr = abs(r2 - r);
        const int dc = abs(c2 - c);

        if (dr > 0 && dc > 0) {
            return fail("Cycle contains a diagonal (cut) move");
        }
        if (dr == 0 && dc == 0) {
            return fail("Cycle contains duplicate consecutive vertex");
        }

        if (dr > 0) {
            if (!(c == left_corridor && c2 == left_corridor) && !(c == right_corridor && c2 == right_corridor)) {
                return fail("Vertical move is not on an external corridor");
            }
        } else { // horizontal
            // For full-row paths, horizontal moves must cross an entire row corridor-to-corridor.
            if (!((c == left_corridor && c2 == right_corridor) || (c == right_corridor && c2 == left_corridor))) {
                return fail("Horizontal move is not a full-row corridor-to-corridor crossing");
            }
            rows_crossed_horizontally.insert(r);
        }
    }

    const int path_cost = compute_cycle_cost(sol.cycle);
    if (path_cost != sol.cost) {
        return fail("Reported cost is inconsistent with the path");
    }

    for (const int r : sol.traversed_rows) {
        if (r < 0 || r >= rows) {
            return fail("Traversed row out of bounds");
        }
        if (rows_crossed_horizontally.find(r) == rows_crossed_horizontally.end()) {
            return fail("A traversed row is not actually crossed in the path");
        }
    }

    // Allow row 0 to be crossed for closure even if not listed among traversed_rows.
    for (const int r : rows_crossed_horizontally) {
        if (r == 0) {
            continue;
        }
        if (find(sol.traversed_rows.begin(), sol.traversed_rows.end(), r) == sol.traversed_rows.end()) {
            return fail("Path crosses a non-listed row with a full-row move");
        }
    }

    return true;
}

bool util::is_partial_row_single_column_solution_feasible(const solution &sol, const int rows, const int internal_cols,
                                                          const int budget, string *reason) {
    auto fail = [&](const string &msg) {
        if (reason != nullptr) {
            *reason = msg;
        }
        return false;
    };

    if (rows <= 0 || internal_cols <= 0) {
        return fail("Invalid instance dimensions");
    }

    const int left_corridor = 0;
    const int right_corridor = internal_cols + 1;
    const int graph_cols = internal_cols + 2;

    if (sol.cost > budget) {
        return fail("Solution cost exceeds budget");
    }

    if (sol.cycle.empty()) {
        if (sol.cost != 0) {
            return fail("Empty cycle with non-zero cost");
        }
        return true;
    }

    if (sol.cycle.front() != make_pair(0, left_corridor) || sol.cycle.back() != make_pair(0, left_corridor)) {
        return fail("Cycle must start/end at depot (0,0)");
    }

    int active_side_corridor = -1;

    for (size_t i = 0; i < sol.cycle.size(); ++i) {
        const auto [r, c] = sol.cycle[i];
        if (r < 0 || r >= rows) {
            return fail("Cycle row out of bounds");
        }
        if (c < 0 || c >= graph_cols) {
            return fail("Cycle column out of bounds");
        }

        if (i + 1 == sol.cycle.size()) {
            continue;
        }

        const auto [r2, c2] = sol.cycle[i + 1];
        const int dr = abs(r2 - r);
        const int dc = abs(c2 - c);

        if (dr == 0 && dc == 0) {
            return fail("Cycle contains duplicate consecutive vertex");
        }
        if (dr > 0 && dc > 0) {
            return fail("Cycle contains a diagonal (cut) move");
        }

        if (dr > 0) {
            if (c != c2) {
                return fail("Vertical move changes column");
            }
            if (c != left_corridor && c != right_corridor) {
                return fail("Vertical move is not on an external corridor");
            }
        } else { // horizontal
            const bool touches_left = (c == left_corridor || c2 == left_corridor);
            const bool touches_right = (c == right_corridor || c2 == right_corridor);
            if (!touches_left && !touches_right) {
                return fail("Horizontal move does not touch an external corridor");
            }

            // Infer the active side corridor from the first non-top-row partial excursion.
            if (r > 0 && active_side_corridor == -1) {
                if (touches_left) {
                    active_side_corridor = left_corridor;
                } else if (touches_right) {
                    active_side_corridor = right_corridor;
                }
            }

            // On rows below the top row, partial-row excursions must use a single side corridor.
            if (r > 0 && active_side_corridor != -1) {
                const bool touches_active = (c == active_side_corridor || c2 == active_side_corridor);
                if (!touches_active) {
                    return fail("Horizontal move uses the wrong side corridor for single-column partial policy");
                }
            }
        }
    }

    const int path_cost = compute_cycle_cost(sol.cycle);
    if (path_cost != sol.cost) {
        return fail("Reported cost is inconsistent with the path");
    }

    return true;
}

bool util::is_partial_row_solution_feasible(const solution &sol, const int rows, const int internal_cols,
                                            const int budget, string *reason) {
    auto fail = [&](const string &msg) {
        if (reason != nullptr) {
            *reason = msg;
        }
        return false;
    };

    if (rows <= 0 || internal_cols <= 0) {
        return fail("Invalid instance dimensions");
    }

    const int left_corridor = 0;
    const int right_corridor = internal_cols + 1;
    const int graph_cols = internal_cols + 2;

    if (sol.cost > budget) {
        return fail("Solution cost exceeds budget");
    }

    if (sol.cycle.empty()) {
        if (sol.cost != 0) {
            return fail("Empty cycle with non-zero cost");
        }
        return true;
    }

    if (sol.cycle.front() != make_pair(0, left_corridor) || sol.cycle.back() != make_pair(0, left_corridor)) {
        return fail("Cycle must start/end at depot (0,0)");
    }

    for (size_t i = 0; i < sol.cycle.size(); ++i) {
        const auto [r, c] = sol.cycle[i];
        if (r < 0 || r >= rows) {
            return fail("Cycle row out of bounds");
        }
        if (c < 0 || c >= graph_cols) {
            return fail("Cycle column out of bounds");
        }

        if (i + 1 == sol.cycle.size()) {
            continue;
        }

        const auto [r2, c2] = sol.cycle[i + 1];
        const int dr = abs(r2 - r);
        const int dc = abs(c2 - c);

        if (dr == 0 && dc == 0) {
            return fail("Cycle contains duplicate consecutive vertex");
        }
        if (dr > 0 && dc > 0) {
            return fail("Cycle contains a diagonal (cut) move");
        }

        if (dr > 0) {
            if (c != c2) {
                return fail("Vertical move changes column");
            }
            if (c != left_corridor && c != right_corridor) {
                return fail("Vertical move is not on an external corridor");
            }
        }
    }

    const int path_cost = compute_cycle_cost(sol.cycle);
    if (path_cost != sol.cost) {
        return fail("Reported cost is inconsistent with the path");
    }

    return true;
}
