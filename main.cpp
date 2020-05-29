#include "matrix.h"
#include "functions.h"

#include <iostream>
#include <iomanip>

using namespace std;

#include <vector>
#include <exception>
#include <fstream>
#include <climits>


int main() {
    Polynom x ({1, 1});

    Matrix<Polynom> a = {{4, -1, 5, 9},
                         {1, 1, 5, 2},
                         {7, -7, 0, -1},
                         {-1, -1, 1, 8}};

    Matrix<Polynom> ones(4, Polynom(1));

    a -= x * ones;

    Polynom det = a.det();

    std::cout << det;

    for (auto i: det.roots()) {
        std::cout << i << std::endl;
    }
}
