#include <iostream>
#include <vector>
#include <cmath>


#include "Lanczos.h"


int main() {
	CSR<double> A = CSR<double>(3, 3);
    A.set(0, 0, 1.0);   A.set(0, 1, 3.0);   A.set(0, 2, 5.0); 
    A.set(1, 0, 3.0);   A.set(1, 1, 7.0);   A.set(1, 2, 10.0); 
    A.set(2, 0, 5.0);   A.set(2, 1, 10.0);  A.set(2, 2, 12.0); 

    CSR<double> B = CSR<double>(3, 3);
    B.set(0, 0, 4.0);   B.set(0, 1, 5.0);   B.set(0, 2, 6.0); 
    B.set(1, 0, 5.0);   B.set(1, 1, 7.0);   B.set(1, 2, 9.0); 
    B.set(2, 0, 6.0);   B.set(2, 1, 9.0);   B.set(2, 2, 5.0); 

	int m = 3;

	std::cout << A << B << std::endl;

	std::vector<double> alpha, beta;
	std::vector<std::vector<double> > q;
	LanczosInversePowerProcessForGeneral(A, B, alpha, beta, q, m);

	for(int i = 0; i < m; i++){
		double lambda = BisectionMethod(alpha, beta, i);
		std::cout << 1.0 / lambda << std::endl;
		
		std::vector<double> y = InversePowerMethod(alpha, beta, lambda);
		std::vector<double> x = ReconvertVector(y, q);
		for(auto xi : x){
			std::cout << xi << "\t";
		}
		std::cout << std::endl;
	}	

	return 0;
}