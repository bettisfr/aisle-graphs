#ifndef DEPLOYMENT_H
#define DEPLOYMENT_H

#include <string>

#include "../io/output.h"

using namespace std;

struct random_instance_params {
    int rows = 10;
    // Number of INTERNAL columns n (profit-bearing only).
    // Corridor columns are implicit and have profit 0.
    int cols = 20;
    int seed = 0;
    int min_reward = 0;
    int max_reward = 10;
};

class deployment {
private:
    string label_;
    RewardGrid rewards_;
    int seed_ = -1;

public:
    deployment() = default;
    deployment(string label, RewardGrid rewards, int seed = -1);

    static deployment from_rewards(string label, const RewardGrid &rewards);
    static deployment random_grid(const random_instance_params &params);

    [[nodiscard]] const string &label() const;
    [[nodiscard]] const RewardGrid &rewards() const;
    [[nodiscard]] int rows() const;
    // Returns n = number of INTERNAL columns stored in RewardGrid.
    [[nodiscard]] int cols() const;
    [[nodiscard]] int seed() const;

    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool is_rectangular() const;
};

#endif //DEPLOYMENT_H
