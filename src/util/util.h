#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <vector>

#include "../io/output.h"

using namespace std;

class util {

public:
    static double compute_cycle_cost(const vector<pair<int, int>> &cycle);
    static vector<pair<int, int>> build_full_row_cycle_from_selection(const vector<int> &selected_rows_0based,
                                                                      int internal_cols);

    static bool is_full_row_solution_feasible(const solution &, int rows, int internal_cols, int budget,
                                              string *reason = nullptr);
    static bool is_partial_row_single_column_solution_feasible(const solution &, int rows, int internal_cols,
                                                               int budget, string *reason = nullptr);

    template<typename T>
    static pair<double, double> calculate_avg_std(const vector<T> &values) {
        if (values.empty()) {
            return {0.0, 0.0};
        }

        T sum = static_cast<T>(0);
        for (const T &value: values) {
            sum += value;
        }

        double average = static_cast<double>(sum) / values.size();

        if (values.size() == 1) {
            return {average, 0.0}; // Standard deviation is 0 for a single value
        }

        T sum_squared_diff = static_cast<T>(0);
        for (const T &value: values) {
            T diff = value - static_cast<T>(average);
            sum_squared_diff += diff * diff;
        }

        double std_dev = sqrt(static_cast<double>(sum_squared_diff) / (values.size() - 1.0));

        return {average, std_dev};
    }
};

#endif //UTIL_H
