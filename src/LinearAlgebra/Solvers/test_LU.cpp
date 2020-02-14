#include <iostream>

#include "../Models/Matrix.h"
#include "LU.h"

using namespace PANSFEM2;

int main(){
    Matrix<double> A = Matrix<double>(4, 4);
    A(0, 0) = 1.0;  A(0, 1) = 1.0;  A(0, 2) = 0.0;  A(0, 3) = 3.0;
    A(1, 0) = 2.0;  A(1, 1) = 1.0;  A(1, 2) = -1.0; A(1, 3) = 1.0;
    A(2, 0) = 3.0;  A(2, 1) = -1.0; A(2, 2) = -1.0; A(2, 3) = 2.0;
    A(3, 0) = -1.0; A(3, 1) = 2.0;  A(3, 2) = 3.0;  A(3, 3) = -1.0;
    std::vector<int> pivot = std::vector<int>(4);
    LU(A, pivot);
    Vector<double> b = Vector<double>(4);
    b(0) = 4.0;
    b(1) = 1.0;
    b(2) = -3.0;
    b(3) = 4.0;
    SolveLU(A, b, pivot);
    std::cout << b << std::endl;
    
    return 0;
}