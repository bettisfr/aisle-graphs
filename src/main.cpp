#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "core/algorithms.h"
#include "core/deployment.h"
#include "io/input.h"
#include "io/output.h"
#include "util/util.h"

using namespace std;
using namespace chrono;

namespace {

deployment build_instance(const input &cfg) {
    if (!cfg.input_csv_path.empty()) {
        RewardGrid grid = load_reward_grid_csv(cfg.input_csv_path);
        return deployment::from_rewards(cfg.input_csv_path, grid);
    }

    random_instance_params params;
    params.rows = cfg.rows;
    params.cols = cfg.cols;
    params.seed = cfg.seed;
    params.min_reward = cfg.min_reward;
    params.max_reward = cfg.max_reward;
    return deployment::random_grid(params);
}

experiment_series run_batch_experiment(algorithms &alg, const deployment &dep, const input &cfg) {
    experiment_series series;
    series.algorithm_key = algorithm_to_string(cfg.algorithm);
    series.instance_label = dep.label();

    if (cfg.budget_points <= 1) {
        solution s = alg.run_experiment(cfg.algorithm, cfg.min_budget);
        series.points.push_back({cfg.min_budget, s.reward, s.cost});
        return series;
    }

    for (int i = 0; i < cfg.budget_points; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(cfg.budget_points - 1);
        const int budget = static_cast<int>(cfg.min_budget + alpha * (cfg.max_budget - cfg.min_budget));
        solution s = alg.run_experiment(cfg.algorithm, budget);
        series.points.push_back({budget, s.reward, s.cost});
    }

    return series;
}

} // namespace

int main(const int argc, char **argv) {
    cout << fixed << setprecision(2);

    try {
        input cfg = read_parameters(argc, argv);
        if (check_parameters(cfg)) {
            return -1;
        }

        if (cfg.log == 1) {
            print_parameters(cfg);
        }

        deployment dep = build_instance(cfg);
        if (dep.empty() || !dep.is_rectangular()) {
            cerr << "Invalid reward grid (empty or non-rectangular)." << endl;
            return -1;
        }

        algorithms alg(dep);

        if (cfg.run_batch == 1) {
            experiment_series series = run_batch_experiment(alg, dep, cfg);
            if (cfg.save == 1) {
                filesystem::create_directories("output");
                const string out_path = "output/" + cfg.exp_name + ".csv";
                save_experiment_series_csv(out_path, series);
                cout << "Saved batch output to " << out_path << endl;
            } else {
                cout << "Batch points generated: " << series.points.size() << endl;
            }
            return 0;
        }

        auto t0 = high_resolution_clock::now();
        solution sol = alg.run_experiment(cfg.algorithm, cfg.budget);
        auto t1 = high_resolution_clock::now();
        sol.running_time_sec = static_cast<double>(duration_cast<microseconds>(t1 - t0).count()) / 1e6;

        if (cfg.algorithm == 0 || cfg.algorithm == 1) {
            string reason;
            const bool ok = util::is_full_row_solution_feasible(sol, dep.rows(), dep.cols(), cfg.budget, &reason);
            cout << "full_row_feasible=" << (ok ? 1 : 0) << "\n";
            cout << "full_row_reason=" << reason << "\n";
        }

        print_solution_summary(sol);
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return -1;
    }

    return 0;
}
