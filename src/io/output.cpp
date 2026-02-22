#include "output.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

using namespace std;

void print_solution_summary(const solution &sol) {
    cout << "algorithm=" << sol.algorithm_key << "\n";
    cout << "reward=" << sol.reward << "\n";
    cout << "cost=" << sol.cost << "\n";
    cout << "runtime_sec=" << sol.running_time_sec << "\n";
    cout << "traversed_rows=" << sol.traversed_rows.size() << "\n";
    cout << "traversed_rows_list=";
    for (size_t i = 0; i < sol.traversed_rows.size(); ++i) {
        if (i > 0) {
            cout << ",";
        }
        cout << sol.traversed_rows[i];
    }
    cout << "\n";
    cout << "cycle_points=" << sol.cycle.size() << "\n";
    cout << "cycle_list=";
    for (size_t i = 0; i < sol.cycle.size(); ++i) {
        if (i > 0) {
            cout << ";";
        }
        cout << sol.cycle[i].first << ":" << sol.cycle[i].second;
    }
    cout << "\n";
    cout << "selected_columns_per_row=" << sol.selected_columns_per_row.size() << "\n";
    if (!sol.notes.empty()) {
        cout << "notes:\n";
        for (const auto &note : sol.notes) {
            cout << "  - " << note << "\n";
        }
    }
}

void save_experiment_series_csv(const string &path, const experiment_series &series) {
    ofstream out(path, ofstream::out | ofstream::trunc);
    if (!out.is_open()) {
        throw runtime_error("Unable to open output file: " + path);
    }

    out << "algorithm,instance,budget,cost,reward\n";
    out << fixed << setprecision(2);
    for (const auto &point : series.points) {
        out << series.algorithm_key << ","
            << series.instance_label << ","
            << point.budget << ","
            << point.cost << ","
            << point.reward << "\n";
    }
}
