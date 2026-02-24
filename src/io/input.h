#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "output.h"

using namespace std;

struct input {
    int log = 1;
    int save = 0;

    // 0 -> ofr
    // 1 -> gfr
    // 2 -> hprgc
    // 3 -> oprsc
    // 4 -> gprsc
    int algorithm = 0;
    int budget = 0;

    // CSV convention: stores only the n internal profit columns.
    // The two corridor columns c_0 and c_{n+1} are implicit and have profit 0.
    string input_csv_path;
    int input_id = -1;

    // Number of rows m and INTERNAL columns n (not the full graph columns n+2).
    int rows = 10;
    int cols = 20;
    int seed = 0;
    int min_reward = 0;
    int max_reward = 10;

    int run_batch = 0;
    int min_budget = 0;
    int max_budget = 0;
    int budget_points = 20;

    string exp_name = "aisle_graphs";
};

input read_parameters(int argc, char **argv);
input load_parameters_from_file(const string &cfg_path, input defaults = {});
bool check_parameters(const input &cfg);
void print_parameters(const input &cfg);
string algorithm_to_string(int algorithm);
int algorithm_from_string(const string &name);

RewardGrid load_reward_grid_csv(const string &path);
void save_reward_grid_csv(const string &path, const RewardGrid &grid);

#endif //INPUT_H
