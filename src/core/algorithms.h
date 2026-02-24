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
    solution make_todo_solution(const string &algorithm_key, int budget, const string &note) const;
    solution opt_partial_row_single_column_right(int budget) const;
    solution opt_partial_row_single_column_right_public(int budget) const; // internal rooted wrapper for right-side OPRSC

public:
    explicit algorithms(const deployment &dep);

    solution run_experiment(int algorithm, int budget);

    solution opt_full_row(int budget) const;                            // ofr
    solution opt_partial_row_single_column(int budget) const;           // oprsc
    solution heuristic_partial_row(int budget) const;                   // hprgc
    solution apx_partial_row(int budget) const;                        // apr
    solution greedy_full_row(int budget) const;                        // gfr
    solution greedy_partial_row_single_column(int budget) const;        // gprsc
};

#endif //ALGORITHMS_H
