#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include "polynom.h"

template<typename Ring=Polynom>
class Matrix {
private:
    std::pair<size_t, size_t> size_;
    std::vector<std::vector <Ring>> data_;
public:
    Matrix() = delete;
    Matrix(size_t, size_t);
    Matrix(const Matrix&);
    Matrix(size_t, Ring elem);
    Matrix(const std::initializer_list<std::initializer_list<Ring>>&);
    ~Matrix();

    Matrix& operator=(const Matrix& right);
    Ring det() const;
    std::pair<size_t, size_t> size() const;

    template<typename Ring_>
    friend const Matrix<Ring_> operator+(const Matrix<Ring_>& left, const Matrix<Ring_>& right);
    template<typename Ring_>
    friend const Matrix<Ring_> operator-(const Matrix<Ring_>& left, const Matrix<Ring_>& right);
    template <unsigned N_, unsigned M_, unsigned K_, typename Ring_>
    friend const Matrix<Ring_> operator*(const Matrix<Ring_>& left, const Matrix<Ring_>& right);///!!!!!!!!!!!!

    template<typename Ring_>
    friend const Matrix<Ring_> operator*(const Matrix<Ring_>& left, const Ring_& right);
    template<typename Ring_>
    friend const Matrix<Ring_> operator*(const Ring_& left, const Matrix<Ring_>& right);

    template<typename Ring_>
    friend Matrix<Ring_>& operator*=(Matrix<Ring_>& left, const Ring_& right);

    template<typename Ring_>
    friend Matrix<Ring_>& operator+=(Matrix<Ring_>& left, const Matrix<Ring_>& right);
    template<typename Ring_>
    friend Matrix<Ring_>& operator-=(Matrix<Ring_>& left, const Matrix<Ring_>& right);
    template <unsigned N_, typename Ring_>
    friend Matrix<Ring_>& operator*=(Matrix<Ring_>& left, const Matrix<Ring_>& right);

    const std::vector<Ring>& operator[](const size_t ind) const;
    std::vector<Ring>& operator[](const size_t ind);
};

///______________________Constructor_Destructor_____________________///
template<typename Ring>
Matrix<Ring>::Matrix(size_t n, size_t m) :
    size_({n, m}) {
    if (n == 0 || m == 0)
        throw std::runtime_error("Bad size for construct (0)");
    data_.resize(n, std::vector<Ring>(m));
}

template<typename Ring>
Matrix<Ring>::Matrix(const Matrix<Ring>& right) :
    size_(right.size_),
    data_(right.data_) {}

template<typename Ring>
Matrix<Ring>::Matrix(size_t n, Ring elem) : Matrix(n, n) {
    for (size_t i = 0; i < n; ++i) {
        data_[i][i] = elem;
    }
}

template<typename Ring>
Matrix<Ring>::Matrix(const std::initializer_list<std::initializer_list<Ring>>& init) : Matrix(init.size(), init.begin()->size()) {
    size_t i = 0,
           j = 0;
    for (auto it = init.begin(); it != init.end(); it++, i++) {
        j = 0;
        if (it->size() != size_.second)
            throw std::runtime_error("Bad initializer list");
        for (auto jt = it->begin(); jt != it->end(); jt++, j++) {
            data_[i][j] = *jt;
        }
    }
}

template<typename Ring>
Matrix<Ring>::~Matrix() {}

///____________________Assignment_____________________________///
template<typename Ring>
Matrix<Ring>& Matrix<Ring>::operator=(const Matrix<Ring>& right) {
    if (this == &right) {
        return *this;
    }
    size_ = right.size_;
    data_ = right.data_;
    return *this;
}

template<typename Ring>
const std::vector<Ring>& Matrix<Ring>::operator[](const size_t ind) const {
    return data_[ind];
}

template<typename Ring>
std::vector<Ring>& Matrix<Ring>::operator[](const size_t ind) {
    return data_[ind];
}

///_______________________IO________________________________________///
template<typename Ring_>
std::istream& operator>>(std::istream &in, Matrix<Ring_> &matrix) {
    for (size_t i = 0; i < matrix.size().first; ++i) {
        for (size_t c = 0; c < matrix.size().second; ++c) {
            in >> matrix[i][c];
        }
    }
    return in;
}

template<typename Ring_>
std::ostream& operator<<(std::ostream &out, const Matrix<Ring_> &matrix) {
    for (size_t i = 0; i < matrix.size().first; ++i) {
        for (size_t c = 0; c < matrix.size().second; ++c) {
            out << matrix[i][c] << ' ';
        }
        out << std::endl;
    }
    return out;
}

///_______________________Binary________________________________________///
template<typename Ring_>
const Matrix<Ring_> operator+(const Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    Matrix<Ring_> sum = left;
    return sum += right;
}

 template<typename Ring_>
const Matrix<Ring_> operator-(const Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    Matrix<Ring_> dif = left;
    return dif -= right;
}

template<typename Ring_>
const Matrix<Ring_> operator*(const Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    auto& l_s = left.size_;
    auto& r_s = right.size_;
    if (l_s.second != r_s.first)
        throw std::runtime_error("Try to multilpicate matrix with bad sizes");
    Matrix<Ring_> multy(l_s.first, r_s.second);
    for (size_t i = 0; i < l_s.first; ++i) {
        for (size_t c = 0; c < r_s.second; ++c) {
            for (size_t k = 0; k < l_s.second; ++k) {
                multy[i][c] += left[i][k] * right[k][c];
            }
        }
    }
    return multy;
}

template<typename Ring_>
const Matrix<Ring_> operator*(const Matrix<Ring_>& left, const Ring_& right) {
    Matrix<Ring_> multy = left;
    return multy *= right;
}

template<typename Ring_>
const Matrix<Ring_> operator*(const Ring_& left, const Matrix<Ring_>& right) {
    Matrix<Ring_> multy = right;
    return multy *= left;
}

///_____________________________________________________________________///
template<typename Ring_>
Matrix<Ring_>& operator+=(Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    if (left.size_ != right.size_)
        throw std::runtime_error("Bad sizes of matrix '+'");
    for (size_t i = 0; i < left.size_.first; ++i) {
        for (size_t c = 0; c < right.size_.second; ++c) {
            left[i][c] += right[i][c];
        }
    }
    return left;
}

template<typename Ring_>
Matrix<Ring_>& operator-=(Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    if (left.size_ != right.size_)
        throw std::runtime_error("Bad sizes of matrix for '+'");
    for (size_t i = 0; i < left.size_.first; ++i) {
        for (size_t c = 0; c < right.size_.second; ++c) {
            left[i][c] -= right[i][c];
        }
    }
    return left;
}

template <typename Ring_>
Matrix<Ring_>& operator*=(Matrix<Ring_>& left, const Matrix<Ring_>& right) {
    if (left.size_ != right.size_)
        throw std::runtime_error("Bad sizes of matrix for '*'");
    left = left * right;
    return left;
}

template<typename Ring_>
Matrix<Ring_>& operator*=(Matrix<Ring_>& left, const Ring_& right) {
    for (size_t i = 0; i < left.size_.first; ++i) {
        for (size_t c = 0; c < left.size_.second; ++c) {
            left[i][c] *= right;
        }
    }
    return left;
}

///_______________________Methods_______________________________________///
template<typename Ring>
Ring Matrix<Ring>::det() const {
    if (size_.first != size_.second)
        throw std::runtime_error("Bad sizes of matrix for det");
    Ring div = 1;
    Matrix<Ring> diag = methodGauss(*this, div);
    Ring det = 1;
    for (size_t i = 0; i < size_.first; ++i) {
        det *= diag[i][i];
    }
    det /= div;
    return det;
}

template<typename Ring>
std::pair<size_t, size_t> Matrix<Ring>::size() const {
    return size_;
}

///_________________________________________________///
template<typename Ring>
void sum_elem(std::vector<Ring>& v1, const std::vector<Ring> v2, const Ring& coef) {
    if (v1.size() != v2.size()) {
        throw std::bad_cast();
    }
    for (size_t i = 0; i < v1.size(); ++i) {
        v1[i] += coef * v2[i];
    }
}

template<typename Ring>
void mult_elem(std::vector<Ring>& v, const Ring& coef) {
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] *= coef;
    }
}

template<typename Ring>
Matrix<Ring> methodGauss(const Matrix<Ring>& matrix, Ring& div) {
    Matrix<Ring> diagonal = matrix;
    for (size_t i = 0; i < std::min(matrix.size().first, matrix.size().second); ++i) {
        bool zero = true;
        for (size_t c = i; c < matrix.size().first; ++c) {
            if (diagonal[c][i] != 0) {
                std::swap(diagonal[i], diagonal[c]);
                if (i != c) div *= -1;
                zero = false;
                break;
            }
        }
        for (size_t c = 0; c < matrix.size().first && !zero; ++c) {
            if (c == i || diagonal[c][i] == 0) {
                continue;
            }
            Ring temp = diagonal[c][i];
            mult_elem(diagonal[c], diagonal[i][i]);
            div *= diagonal[i][i];
            sum_elem(diagonal[c], diagonal[i], -temp);
        }
    }
    return diagonal;
}
#endif // MATRIX_H
