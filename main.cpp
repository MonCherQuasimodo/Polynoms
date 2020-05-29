#include "matrix.h"
#include "functions.h"

#include <iostream>

using namespace std;


int main() {
    Polynom x ({1, 1});

    Matrix<Polynom> a = {{3, -8, 16, 2},
                         {7, 1, 95, 0},
                         {7, -7, 0, -1},
                         {-1, -1, 1, 8}};

    ///______________Можно добавить ввод из файла!!!_____________///
    Matrix<Polynom> ones(4, Polynom(1));

    std::cout << "Изначальная матрица" << std::endl;
    std::cout << a << std::endl;
    a -= x * ones;
    std::cout << "A - x * E" << std::endl;
    std::cout << a << std::endl;
    Polynom det = a.det();
    std::cout << "Характеристический многочлен" << std::endl;
    std::cout << det << std::endl << std::endl;
    std::cout << "Собственные значения:" << std::endl;

    for (auto i: det.roots()) {
        std::cout << i << std::endl;
    }

}
