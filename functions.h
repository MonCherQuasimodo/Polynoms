#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdexcept>
#include <deque>
#include <set>
#include <type_traits>
#include <algorithm>
#include <iterator>

const double EPSILON = 1e-10;

template <typename Function>
double newton_method(const Function& f, const Function& f_diff, const Function& f_diff_2, double x_l, double x_r) {
    if (f(x_l) == 0 || f(x_r) == 0)
        return f(x_l) == 0 ? x_l : x_r;
    if (f(x_l) * f(x_r) > 0)
        throw std::runtime_error("Bad interval");

    double x_med = (x_l + x_r) / 2;
    double x_0 = f_diff(x_med) * f_diff_2(x_med) > 0 ? x_r : x_l;
    double x_1 = x_0 - f(x_0) / f_diff(x_0);
    while (std::abs(x_0 - x_1) > EPSILON) {
        x_0 = x_1;
        x_1 = x_1 - f(x_1) / f_diff(x_1);
    }
    return x_1;
}

enum SIDE {
    LEFT = -1,
    RIGHT = 1
};

template <typename Function>
double find_border_point_from_last_extremum_point(const Function& f, double x, SIDE side) {
    double x_add = 2 * static_cast<int>(side);
    if (f.differenced()(x + x_add) * f(x) * static_cast<int>(side) > 0 || f(x) == 0)
        throw std::runtime_error("Bad point");
    const double value = f(x);
    while (value * f(x + x_add) > 0 ) {
        x_add *= 2;
    }
    return x + x_add;
}

template <typename Function>
std::set<double> find_roots(const Function& func) {
    std::deque<Function> all_diff;
    auto func_copy = func;
    for (size_t i = 0; i <= static_cast<size_t>(func.degree()); ++i) {
        all_diff.push_front(func_copy);
        func_copy.difference();
    }

    std::deque<std::set<double>> zeros(all_diff.size());
    zeros[1] = all_diff[1].roots();
    for (size_t i = 2; i <= static_cast<size_t>(func.degree()); ++i) {
        std::set<double> stab_points;
        stab_points.insert(zeros[i-1].begin(), zeros[i-1].end());
        stab_points.insert(zeros[i-2].begin(), zeros[i-2].end());
        if (stab_points.empty()) {
            stab_points.insert(0);
        }
        try {
            stab_points.insert(find_border_point_from_last_extremum_point(all_diff[i],
                                                                          *stab_points.begin(), LEFT));
        }  catch (...) {};
        try {
            stab_points.insert(find_border_point_from_last_extremum_point(all_diff[i],
                                                                          *stab_points.rbegin(), RIGHT));
        }  catch (...) {};
        for (auto it_l = stab_points.begin(), it_r = std::next(it_l); it_r != stab_points.end(); ++it_l, ++it_r) {
            if (all_diff[i](*it_l) * all_diff[i](*it_r) < 0) {
                zeros[i].insert(newton_method(all_diff[i], all_diff[i-1], all_diff[i-2], *it_l, *it_r));
            }
            if (all_diff[i](*it_l) == 0) {
                zeros[i].insert(*it_l);
            }
            if (all_diff[i](*it_r) == 0) {
                zeros[i].insert(*it_r);
                ++it_l;
                ++it_r;
            }
        }
    }
    return zeros.back();
}
#endif // FUNCTIONS_H
