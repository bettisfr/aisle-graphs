#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <iostream>

using namespace std;

class util {

public:
    static double get_2D_distance(double, double, double, double);

    static double get_3D_distance(double, double, double, double, double, double);

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
