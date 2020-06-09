#pragma once

#include "functions.h"

#include <iostream>
#include <list>
#include <set>
#include <cmath>
#include <stdexcept>

class Polynom {
private:
    using Monom = std::pair<long double, size_t>;
    std::list<Monom> data_;
    static Polynom& alg_sum_(Polynom& left, const Polynom& right, bool substract=false);
    static std::pair<Polynom, Polynom> division_(const Polynom& left, const Polynom& right);
public:
    Polynom(const long double&);
    Polynom(const Monom&);

    Polynom() = default;
    Polynom(const Polynom&) = default;
    Polynom(Polynom&&) = default;
    Polynom& operator=(const Polynom&) = default;
    Polynom& operator=(Polynom&&) = default;
    ~Polynom() = default;

    explicit operator int() const;
    explicit operator long double() const;
    long double operator()(const long double&) const;

    int degree() const;
    void clear();

    std::set<long double> roots(METHOD=NEWTONE) const;

    void difference(size_t n=1);
    const Polynom differenced(size_t n=1) const;

    const Polynom operator-() const;

    friend Polynom& operator+=(Polynom& left, const Polynom& right);
    friend Polynom& operator-=(Polynom& left, const Polynom& right);
    friend Polynom& operator*=(Polynom& left, const Polynom& right);
    friend Polynom& operator/=(Polynom& left, const Polynom& right);
    friend Polynom& operator%=(Polynom& left, const Polynom& right);

    friend const Polynom operator+(const Polynom& left, const Polynom& right);
    friend const Polynom operator-(const Polynom& left, const Polynom& right);
    friend const Polynom operator*(const Polynom& left, const Polynom& right);
    friend const Polynom operator/(const Polynom& left, const Polynom& right);
    friend const Polynom operator%(const Polynom& left, const Polynom& right);

    friend std::istream& operator>>(std::istream& in, Polynom& number);
    friend std::ostream& operator<<(std::ostream& out, const Polynom& number);

    friend bool operator==(const Polynom& left, const Polynom& right);
    friend bool operator!=(const Polynom& left, const Polynom& right);
};


Polynom::Polynom(const long double& scalar) :
    Polynom({scalar, 0}) {}

Polynom::Polynom(const Monom& monom) {
    if (monom.first != 0) {
        data_.push_back(monom);
    }
}

Polynom::operator long double() const {
    if (data_.size() == 1 && data_.front().second == 0) {
        return data_.front().first;
    }
    throw std::bad_cast();
}

Polynom::operator int() const {
    long double result = static_cast<long double>(*this);
    if (result == static_cast<int>(result)) {
        return result;
    }
    throw std::bad_cast();
}

long double Polynom::operator()(const long double& x) const {
    long double result = 0;
    for (auto i: data_) {
        result += i.first * std::pow(x, i.second);
    }
    return result;
}

int Polynom::degree() const {
    return data_.size() ? data_.back().second : -1;
}

void Polynom::clear() {
    data_.clear();
}

std::set<long double> Polynom::roots(METHOD method) const {
    switch (degree()) {
    case -1:
    case 0:
        return {};
    case 1: {
        auto& d = data_;
        return {d.size() == 1 ? 0 : -d.front().first / d.back().first};
    }
    default:
        return find_roots(*this, method);
    }
}

void Polynom::difference(size_t n) {
    for (size_t c = 0; c < n; ++c) {
        if (data_.size() && data_.front().second == 0) {
            data_.pop_front();
        }
        for (auto& i: data_) {
            i.first *= i.second;
            --i.second;
        }
    }
}

const Polynom Polynom::differenced(size_t n) const {
    Polynom result = (*this);
    result.difference(n);
    return result;
}

const Polynom Polynom::operator-() const {
    return (-1) * (*this);
}

Polynom& Polynom::alg_sum_(Polynom& left, const Polynom& right, bool substract) {
    auto& l = left.data_;
    auto& r = right.data_;
    auto it_l = l.begin();
    auto it_r = r.begin();
    long double coef = (substract ? -1 : 1);
    while (it_l != l.end() && it_r != r.end()) {
        if (it_l->second > it_r->second) {
            l.insert(it_l, {it_r->first * coef, it_r->second});
            ++it_r;
        } else if (it_l->second == it_r->second) {
            it_l->first += coef * it_r->first;
            it_l = (it_l->first == 0 ? l.erase(it_l) : ++it_l);
            ++it_r;
        } else {
            ++it_l;
        }
    }
    while (it_r != r.end()) {
        l.push_back({it_r->first * coef, it_r->second});
        ++it_r;
    }
    return left;
}

std::pair<Polynom, Polynom> Polynom::division_(const Polynom& left, const Polynom& right) {
    if (right.degree() == -1) {
        throw std::runtime_error("Division by 0");
    }

    Polynom mod = left;
    Polynom div = 0;

    auto& m = mod.data_;
    auto& r = right.data_;
    auto& d = div.data_;

    while (mod.degree() >= right.degree()) {
        auto& m_b = m.back();
        auto& r_b = r.back();
        d.push_front({m_b.first / r_b.first, m_b.second - r_b.second});
        mod -= right * d.front();
    }
    return {div, mod};
}

const Polynom operator+(const Polynom& left, const Polynom& right) {
    Polynom result = left;
    return result += right;
}

const Polynom operator-(const Polynom& left, const Polynom& right) {
    Polynom result = left;
    return result -= right;
}

const Polynom operator*(const Polynom& left, const Polynom& right) {
    Polynom result = left;
    return result *= right;
}

const Polynom operator/(const Polynom& left, const Polynom& right) {
    Polynom result = left;
    return result /= right;
}

const Polynom operator%(const Polynom& left, const Polynom& right) {
    Polynom result = left;
    return result %= right;
}


Polynom& operator+=(Polynom& left, const Polynom& right) {
    return Polynom::alg_sum_(left, right);
}

Polynom& operator-=(Polynom& left, const Polynom& right) {
    return Polynom::alg_sum_(left, right, true);
}

Polynom& operator*=(Polynom& left, const Polynom& right){
    Polynom result;
    for (auto i: left.data_) {
        Polynom temp;
        for (auto j: right.data_) {
            temp.data_.push_back({i.first * j.first, i.second + j.second});
        }
        result += temp;
    }
    return left = result;
}

Polynom& operator/=(Polynom& left, const Polynom& right){
    return left = Polynom::division_(left, right).first;
}

Polynom& operator%=(Polynom& left, const Polynom& right){
    return left = Polynom::division_(left, right).second;
}

std::istream& operator>>(std::istream& in, Polynom& poly) {
    //Упрощение для данной задачи.
    poly.clear();
    long double c;
    in >> c;
    poly = Polynom(c);
    return in;
}

std::ostream& operator<<(std::ostream& out, const Polynom& number) {
    auto& data = number.data_;
    out << "{";
    if (data.size() == 0) {
        out << "0";
    }
    for (auto i = data.rbegin(); i != data.rend(); ++i) {
        if (i != data.rbegin()) {
            out << (i->first > 0 ? '+' : '-') << ' ';
        } else if (i->first < 0) {
            out << '-';
        }
        out << std::abs(i->first) << "x^" << i->second << ' ';
    }
    out << "}";
    return out;
}

bool operator==(const Polynom& left, const Polynom& right) {
    return left.data_ == right.data_;
}

bool operator!=(const Polynom& left, const Polynom& right) {
    return !(left == right);
}

const Polynom pow(const Polynom& base, const Polynom& x) {
    size_t n = static_cast<int>(x);
    Polynom a = base;
    Polynom result = 1;
    while (n) {
        if (n & 1)
            result *= a;
        a *= a;
        n >>= 1;
    }
    return result;
}
