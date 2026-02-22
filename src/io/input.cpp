#include "input.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace {
bool is_flag(const string &arg) {
    return !arg.empty() && arg[0] == '-';
}
}

string algorithm_to_string(const int algorithm) {
    switch (algorithm) {
        case 0: return "ofr";
        case 1: return "oprsc";
        case 2: return "oprsc-nv";
        case 3: return "hprgc";
        case 4: return "gfr";
        case 5: return "gprsc";
        default: return "unknown";
    }
}

int algorithm_from_string(const string &name) {
    if (name == "ofr") return 0;
    if (name == "oprsc") return 1;
    if (name == "oprsc-nv") return 2;
    if (name == "hprgc") return 3;
    if (name == "gfr") return 4;
    if (name == "gprsc") return 5;
    throw invalid_argument("Unknown algorithm name: " + name);
}

input load_parameters_from_file(const string &cfg_path, input cfg) {
    ifstream in(cfg_path);
    if (!in.is_open()) {
        throw runtime_error("Unable to open config file: " + cfg_path);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        const size_t pos = line.find('=');
        if (pos == string::npos) {
            continue;
        }

        const string key = line.substr(0, pos);
        const string value = line.substr(pos + 1);

        if (key == "log") cfg.log = stoi(value);
        else if (key == "save") cfg.save = stoi(value);
        else if (key == "algorithm") {
            if (!value.empty() && isdigit(static_cast<unsigned char>(value[0]))) {
                cfg.algorithm = stoi(value);
            } else {
                cfg.algorithm = algorithm_from_string(value);
            }
        }
        else if (key == "budget") cfg.budget = stoi(value);
        else if (key == "input_csv_path") cfg.input_csv_path = value;
        else if (key == "input_id") cfg.input_id = stoi(value);
        else if (key == "rows") cfg.rows = stoi(value);
        else if (key == "cols") cfg.cols = stoi(value);
        else if (key == "seed") cfg.seed = stoi(value);
        else if (key == "min_reward") cfg.min_reward = stoi(value);
        else if (key == "max_reward") cfg.max_reward = stoi(value);
        else if (key == "run_batch") cfg.run_batch = stoi(value);
        else if (key == "min_budget") cfg.min_budget = stoi(value);
        else if (key == "max_budget") cfg.max_budget = stoi(value);
        else if (key == "budget_points") cfg.budget_points = stoi(value);
        else if (key == "exp_name") cfg.exp_name = value;
    }

    return cfg;
}

input read_parameters(const int argc, char **argv) {
    input cfg;

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];

        if (arg == "--file") {
            if (i + 1 >= argc) {
                throw invalid_argument("Missing cfg path after --file");
            }
            cfg = load_parameters_from_file(argv[++i], cfg);
            continue;
        }

        if (!is_flag(arg)) {
            continue;
        }
        if (i + 1 >= argc) {
            throw invalid_argument("Missing value for flag " + arg);
        }

        const string value = argv[++i];

        if (arg == "-algorithm") {
            if (!value.empty() && isdigit(static_cast<unsigned char>(value[0]))) {
                cfg.algorithm = stoi(value);
            } else {
                cfg.algorithm = algorithm_from_string(value);
            }
        }
        else if (arg == "-budget") cfg.budget = stoi(value);
        else if (arg == "-input_csv") cfg.input_csv_path = value;
        else if (arg == "-input_id") cfg.input_id = stoi(value);
        else if (arg == "-rows") cfg.rows = stoi(value);
        else if (arg == "-cols") cfg.cols = stoi(value);
        else if (arg == "-seed") cfg.seed = stoi(value);
        else if (arg == "-min_reward") cfg.min_reward = stoi(value);
        else if (arg == "-max_reward") cfg.max_reward = stoi(value);
        else if (arg == "-run_batch") cfg.run_batch = stoi(value);
        else if (arg == "-min_budget") cfg.min_budget = stoi(value);
        else if (arg == "-max_budget") cfg.max_budget = stoi(value);
        else if (arg == "-budget_points") cfg.budget_points = stoi(value);
        else if (arg == "-save") cfg.save = stoi(value);
        else if (arg == "-log") cfg.log = stoi(value);
        else if (arg == "-exp_name") cfg.exp_name = value;
        else {
            throw invalid_argument("Unknown option: " + arg);
        }
    }

    return cfg;
}

bool check_parameters(const input &cfg) {
    bool error = false;

    if (cfg.rows <= 0) cerr << "Error: rows must be > 0" << endl, error = true;
    if (cfg.cols <= 0) cerr << "Error: cols must be > 0" << endl, error = true;
    if (cfg.algorithm < 0 || cfg.algorithm > 5) cerr << "Error: algorithm must be in [0, 5]" << endl, error = true;
    if (cfg.min_reward > cfg.max_reward) cerr << "Error: min_reward > max_reward" << endl, error = true;
    if (cfg.budget < 0) cerr << "Error: budget must be >= 0" << endl, error = true;
    if (cfg.budget_points <= 0) cerr << "Error: budget_points must be > 0" << endl, error = true;
    if (cfg.run_batch && cfg.max_budget < cfg.min_budget) cerr << "Error: max_budget < min_budget" << endl, error = true;

    return error;
}

void print_parameters(const input &cfg) {
    cout << "algorithm=" << cfg.algorithm << " (" << algorithm_to_string(cfg.algorithm) << ")" << endl;
    cout << "budget=" << cfg.budget << endl;
    cout << "input_csv_path=" << (cfg.input_csv_path.empty() ? "<random>" : cfg.input_csv_path) << endl;
    cout << "input_id=" << cfg.input_id << endl;
    cout << "rows=" << cfg.rows << endl;
    cout << "cols=" << cfg.cols
         << " (internal columns only; graph has n+2 columns with zero-profit corridors)"
         << endl;
    cout << "seed=" << cfg.seed << endl;
    cout << "min_reward=" << cfg.min_reward << endl;
    cout << "max_reward=" << cfg.max_reward << endl;
    cout << "run_batch=" << cfg.run_batch << endl;
    cout << "min_budget=" << cfg.min_budget << endl;
    cout << "max_budget=" << cfg.max_budget << endl;
    cout << "budget_points=" << cfg.budget_points << endl;
    cout << "save=" << cfg.save << endl;
    cout << "exp_name=" << cfg.exp_name << endl;
}

RewardGrid load_reward_grid_csv(const string &path) {
    ifstream in(path);
    if (!in.is_open()) {
        throw runtime_error("Unable to open reward grid CSV: " + path);
    }

    RewardGrid grid; // m x n matrix of internal-column profits (corridor columns omitted)
    string line;
    while (getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        vector<double> row;
        string cell;
        stringstream ss(line);
        while (getline(ss, cell, ',')) {
            if (!cell.empty()) {
                row.push_back(stod(cell));
            }
        }

        if (!row.empty()) {
            grid.push_back(move(row));
        }
    }

    return grid;
}

void save_reward_grid_csv(const string &path, const RewardGrid &grid) {
    ofstream out(path);
    if (!out.is_open()) {
        throw runtime_error("Unable to write reward grid CSV: " + path);
    }

    out << fixed << setprecision(2);
    for (const auto &row : grid) {
        for (size_t j = 0; j < row.size(); ++j) {
            if (j > 0) {
                out << ",";
            }
            out << row[j];
        }
        out << "\n";
    }
}
