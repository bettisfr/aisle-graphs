#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <functional>
#include <vector>

#include "deployment.h"
#include "../io/output.h"

using namespace std;

class algorithms {
private:
    const deployment &dep;
    vector<function<solution(algorithms &, int)>> algorithm_functions;
    solution solve_opr_ilp_gurobi(int budget) const; // internal ILP solver entrypoint (TODO)
    solution solve_oprsc_left_on_grid(const RewardGrid &rewards, int budget, const string &algorithm_key,
                                      bool disable_first_row_partial) const;
    solution opt_partial_row_single_column_right(int budget) const;
    solution opt_partial_row_single_column_right_public(int budget) const; // internal rooted wrapper for right-side OPRSC
    solution heuristic_partial_row_impl(int budget, int component_mode) const; // 0=all,1=S1,2=S2,3=S3
    solution heuristic_partial_row_s1(int budget) const;                // internal debug helper (S1)
    solution heuristic_partial_row_s2(int budget) const;                // internal debug helper (S2)
    solution heuristic_partial_row_s3(int budget) const;                // internal debug helper (S3)

public:
    explicit algorithms(const deployment &dep);

    solution run_experiment(int algorithm, int budget);

    // full row
    solution opt_full_row(int budget) const;                            // ofr
    solution greedy_full_row(int budget) const;                         // gfr

    // partial row, single column
    solution opt_partial_row_single_column(int budget) const;           // oprsc
    solution greedy_partial_row_single_column(int budget) const;        // gprsc

    // partial row
    solution opt_partial_row(int budget) const;                         // opr
    solution opt_partial_row_dp(int budget) const;                      // opr-dp
    solution apx_partial_row(int budget) const;                         // apr
    solution heuristic_partial_row(int budget) const;                   // hpr
    solution greedy_partial_row(int budget) const;                      // gpr

};

#endif //ALGORITHMS_H
