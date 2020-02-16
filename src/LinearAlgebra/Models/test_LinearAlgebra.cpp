#include <iostream>
#include "Vector.h"
#include "Matrix.h"


using namespace PANSFEM2;


int main(){
    Matrix<double> A = Matrix<double>(3, 3);
    A(0, 0) = 1.0;  A(0, 1) = 2.0;  A(0, 2) = 3.0;
    A(1, 0) = 4.0;  A(1, 1) = 5.0;  A(1, 2) = 6.0;
    A(2, 0) = 7.0;  A(2, 1) = 8.0;  A(2, 2) = 9.0;

    Vector<double> b = A.Block(2, 0, 1, 3).Transpose();

    std::cout << b << std::endl;

    Vector<double> a = { 1.0, 2.0, 0.0 };
    Vector<double> c = { 0.0, 1.0, -1.0 };

    std::cout << VectorProduct(a, c) << std::endl;

    return 0;
}