#include <iostream>

#include "../Models/CSR.h"
#include "CG.h"

int main(){
    CSR<double> A = CSR<double>(4, 4);
    A.set(0, 0, 1.0);  A.set(0, 1, 1.0);  A.set(0, 2, 0.0);  A.set(0, 3, 3.0);
    A.set(1, 0, 2.0);  A.set(1, 1, 1.0);  A.set(1, 2, -1.0); A.set(1, 3, 1.0);
    A.set(2, 0, 3.0);  A.set(2, 1, -1.0); A.set(2, 2, -1.0); A.set(2, 3, 2.0);
    A.set(3, 0, -1.0); A.set(3, 1, 2.0);  A.set(3, 2, 3.0);  A.set(3, 3, -1.0);
    std::vector<double> x = std::vector<double>(4, 1.0);
    std::vector<double> b = A*x;

    x = BiCGSTAB(A, b, 1000, 1.0e-1);

    for(auto xi : x) {
        std::cout << xi << std::endl;
    } 
    
    return 0;
}