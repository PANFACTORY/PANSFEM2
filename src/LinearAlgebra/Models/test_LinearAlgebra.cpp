#include <iostream>
#include "Vector.h"
#include "Matrix.h"


using namespace PANSFEM2;


int main(){
    Matrix<double> A = Matrix<double>(3, 2);
    A(0, 0) = 1.0;  A(0, 1) = 2.0;
    A(1, 0) = 3.0;  A(1, 1) = 4.0;
    A(2, 0) = 5.0;  A(2, 1) = 6.0;

    Vector<double> b = { 7.0, 8.0 };

    Vector<double> c = A*b;

    std::cout << c << std::endl;

    return 0;
}