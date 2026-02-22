#include <cmath>
#include <algorithm>

#include "util.h"

double util::get_2D_distance(double x1, double y1, double x2, double y2) {
    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

double util::get_3D_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double delta_x = x2 - x1;
    double delta_y = y2 - y1;
    double delta_z = z2 - z1;

    return sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}
