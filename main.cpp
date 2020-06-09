#include "matrix.h"
#include "functions.h"

#include <iostream>

using namespace std;


int main() {
    /*Polynom x ({1, 1});

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
    }*/

    Polynom x ({1, 1});
    Polynom q = 5 * pow(x, 7) + 8 * pow(x, 6) - 78 * pow(x,5) - 78 * pow(x,4) +33 * pow(x, 3) - pow(x,2) + x + 1;
    //Polynom t = 210 * pow(x, 5) + 240 * pow(x, 4) - 1560 * pow(x, 3) - 936 * pow(x, 2) + 198 * x - 2;

    //std::cout << combined_method(t, t.differenced(1), t.differenced(2), -2, 0) << std::endl;


    std::cout << q << std::endl;
    for (auto i: q.roots(NEWTONE)) {
        std::cout << i << std::endl;
    }
    std::cout << std::endl;
    for (auto i: q.roots(SECANT)) {
        std::cout << i << std::endl;
    }
    std::cout << std::endl;
    for (auto i: q.roots(COMBINED)) {
        std::cout << i << std::endl;
    }


}
