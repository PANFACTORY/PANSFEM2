#include <iostream>

#include "../Models/Matrix.h"
#include "LU.h"

using namespace PANSFEM2;

int main(){
    Matrix<double> A = Matrix<double>(3, 3);
    A(0, 0) = 2.0;  A(0, 1) = 3.0;  A(0, 2) = -1.0;
    A(1, 0) = 4.0;  A(1, 1) = 4.0;  A(1, 2) = -3.0;
    A(2, 0) = -2.0; A(2, 1) = 3.0;  A(2, 2) = -1.0;
    LU(A);
    std::cout << A << std::endl;

    Matrix<double> B = Matrix<double>(3, 3);
    B(0, 0) = 3.0;  B(0, 1) = 1.0;  B(0, 2) = 0.0;
    B(1, 0) = 6.0;  B(1, 1) = 1.0;  B(1, 2) = -2.0;
    B(2, 0) = -3.0; B(2, 1) = 0.0;  B(2, 2) = 3.0;
    LU(B);
    std::cout << B << std::endl;

    return 0;
}