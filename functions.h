#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdexcept>
#include <deque>
#include <set>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <climits>

const long double EPSILON = 1e-10;
const long double THE_SECOND_APPROXIMATION = 1.75; ///Const between 1 and 2

template <typename Function>
long double newton_method(const Function& f, const Function& f_diff, const Function& f_diff_2, long double x_l, long double x_r) {
    if (f(x_l) == 0 || f(x_r) == 0)
        return f(x_l) == 0 ? x_l : x_r;
    if (f(x_l) * f(x_r) > 0)
        throw std::runtime_error("Bad interval");
    long double x_med = (x_l + x_r) / 2;
    long double x_0 = f_diff(x_med) * f_diff_2(x_med) > 0 ? x_r : x_l;
    long double x_1 = x_0 - f(x_0) / f_diff(x_0);
    while (std::abs(x_0 - x_1) > EPSILON) {
        x_0 = x_1;
        x_1 = x_1 - f(x_1) / f_diff(x_1);
    }
    return x_1;
}

template <typename Function>
long double secant_method(const Function& f, long double x_1, long double x_2) {
    if (f(x_1) == 0 || f(x_2) == 0)
        return f(x_1) == 0 ? x_1 : x_2;
    if (f(x_1) * f(x_2) < 0)
        return std::abs(x_1 - x_2) < 2 * EPSILON ? (x_1 + x_2) / 2  : throw std::runtime_error("Bad interval");
    while (std::abs(x_1 - x_2) > EPSILON) {
        long double temp = x_2;
        x_2 = x_1 - f(x_1) * (x_2 - x_1) / (f(x_2) - f(x_1));
        x_1 = temp;
    }
    return x_2;
}

template <typename Function>
long double combined_method(const Function& f, const Function& f_diff, const Function& f_diff_2, long double x_l, long double x_r) {
    if (f(x_l) == 0 || f(x_r) == 0)
        return f(x_l) == 0 ? x_l : x_r;
    if (f(x_l) * f(x_r) > 0)
        throw std::runtime_error("Bad interval");
    long double x_med = (x_l + x_r) / 2;
    long double x_b = f_diff(x_med) * f_diff_2(x_med) > 0 ? x_r : x_l;
    long double x_a = x_b == x_r ? x_l : x_r;
    while (std::abs(x_a - x_b) > EPSILON) {
        x_a = x_a - f(x_a) * (x_b - x_a) / (f(x_b) - f(x_a));
        x_b = x_b - f(x_b) / f_diff(x_b);
    }
    return x_a;
}

enum SIDE {
    LEFT = -1,
    RIGHT = 1
};

enum METHOD {
    NEWTONE,
    SECANT,
    COMBINED
};

template <typename Function>
long double find_border_point_from_last_extremum_point(const Function& f, long double x, SIDE side) {
    long double x_add = 2 * static_cast<int>(side);
    if (f.differenced()(x + x_add) * f(x) * static_cast<int>(side) > 0 || f(x) == 0)
        throw std::runtime_error("Bad point");
    const long double value = f(x);
    while (value * f(x + x_add) > 0 ) {
        x_add *= 2;
    }
    return x + x_add;
}

#include <iostream>

template <typename Function>
std::set<long double> find_roots(const Function& func, METHOD method) {
    std::deque<Function> all_diff;
    auto func_copy = func;
    for (size_t i = 0; i <= static_cast<size_t>(func.degree()); ++i) {
        all_diff.push_front(func_copy);
        func_copy.difference();
    }

    std::deque<std::set<long double>> zeros(all_diff.size());
    zeros[1] = all_diff[1].roots();
    for (size_t i = 2; i <= static_cast<size_t>(func.degree()); ++i) {
        std::set<long double> stab_points;
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
                long double root = 0;
                switch (method) {
                case NEWTONE:
                    root = newton_method(all_diff[i], all_diff[i-1], all_diff[i-2], *it_l, *it_r);
                    break;
                case COMBINED:
                    root = combined_method(all_diff[i], all_diff[i-1], all_diff[i-2], *it_l, *it_r);
                    break;
                case SECANT:
                    SIDE side = std::abs(all_diff[i-1](*it_l)) - std::abs(all_diff[i-1](*it_r))  < 0 ? RIGHT : LEFT;
                    long double x_1 = side == RIGHT ? *it_r : *it_l;
                    long double x_2 = x_1 - THE_SECOND_APPROXIMATION * static_cast<int>(side) * EPSILON;
                    root = secant_method(all_diff[i], x_1, x_2);
                    break;
                }
                zeros[i].insert(root);
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
