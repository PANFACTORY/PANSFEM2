#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/FEM/Equation/Numeric.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"
#include "../../src/FEM/Controller/NewtonCotesIntegration.h"


using namespace PANSFEM2;


double f(Vector<double> _x){
    return 4.0/(1.0 + pow(_x(0), 2.0));
}


int main() {
    int n = 50;

    //----------Add nodes----------
	std::vector<Vector<double> > x = std::vector<Vector<double> >(n);
    for(int i = 0; i < n; i++){
        x[i] = { i/(double)(n - 1) };
    }

    //----------Add elements----------
    std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(n - 1);
    for(int i = 0; i < n - 1; i++){
        elements[i] = { i, i + 1 };
    }

    //----------Culculate element value----------
    double value = 0.0;
    for(auto element : elements){
       value += NumericIntegrationOnLine<double, ShapeFunction2Line, NewtonCotes3Line, double(Vector<double>)>(x, element, f);
    }

    std::cout << value << std::endl;
	return 0;
}