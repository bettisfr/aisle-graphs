#include "deployment.h"

#include <random>
#include <stdexcept>
#include <utility>

using namespace std;

deployment::deployment(string label, RewardGrid rewards, int seed)
    : label_(move(label)), rewards_(move(rewards)), seed_(seed) {}

deployment deployment::from_rewards(string label, const RewardGrid &rewards) {
    return deployment(move(label), rewards, -1);
}

deployment deployment::random_grid(const random_instance_params &params) {
    if (params.rows <= 0 || params.cols <= 0) {
        throw invalid_argument("rows and cols must be > 0");
    }
    if (params.min_reward > params.max_reward) {
        throw invalid_argument("min_reward cannot be > max_reward");
    }

    mt19937 gen(params.seed);
    uniform_int_distribution<int> dist(params.min_reward, params.max_reward);

    RewardGrid rewards(params.rows, vector<double>(params.cols, 0.0));
    for (int i = 0; i < params.rows; ++i) {
        for (int j = 0; j < params.cols; ++j) {
            rewards[i][j] = static_cast<double>(dist(gen));
        }
    }

    return deployment("random", move(rewards), params.seed);
}

const string &deployment::label() const { return label_; }
const RewardGrid &deployment::rewards() const { return rewards_; }
int deployment::rows() const { return static_cast<int>(rewards_.size()); }
int deployment::cols() const { return rewards_.empty() ? 0 : static_cast<int>(rewards_.front().size()); }
int deployment::seed() const { return seed_; }
bool deployment::empty() const { return rewards_.empty() || rewards_.front().empty(); }

bool deployment::is_rectangular() const {
    if (rewards_.empty()) {
        return true;
    }
    const size_t width = rewards_.front().size();
    for (const auto &row : rewards_) {
        if (row.size() != width) {
            return false;
        }
    }
    return true;
}
