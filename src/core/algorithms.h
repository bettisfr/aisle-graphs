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

public:
    explicit algorithms(const deployment &dep);

    solution run_experiment(int algorithm, int budget);

    solution opt_full_row(int budget) const;                            // ofr
    solution opt_partial_row_single_column(int budget) const;           // oprsc
    solution heuristic_partial_row(int budget) const;                   // hprgc
    solution greedy_full_row(int budget) const;                        // gfr
    solution greedy_partial_row_single_column(int budget) const;        // gprsc
};

#endif //ALGORITHMS_H
