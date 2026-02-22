#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using RewardGrid = std::vector<std::vector<double>>;

struct solution {
    std::string algorithm_key;

    double reward = 0.0;
    double cost = 0.0;
    double running_time_sec = 0.0;

    std::vector<std::pair<int, int>> cycle;
    std::vector<int> traversed_rows;
    std::vector<int> selected_columns_per_row;
    std::vector<std::string> notes;
};

struct budget_result {
    int budget = 0;
    double reward = 0.0;
    double cost = 0.0;
};

struct experiment_series {
    std::string algorithm_key;
    std::string instance_label;
    std::vector<budget_result> points;
};

void print_solution_summary(const solution &sol);
void save_experiment_series_csv(const std::string &path, const experiment_series &series);

#endif //OUTPUT_H
