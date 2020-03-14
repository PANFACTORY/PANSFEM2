#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/NewtonCotesIntegration.h"


using namespace PANSFEM2;


int main() {
	double f = 0.0;
    double sigma = 0.01;
    double mu = 0.0;

    for(int i = 0; i < NewtonCotes3Line<double>::N; i++){
        Vector<double> x = NewtonCotes3Line<double>::Points[i];
        std::vector<double> w = NewtonCotes3Line<double>::Weights[i];
        //f += w[0]*exp(-0.5*pow((x(0) - mu)/sigma, 2.0))/(sqrt(2.0*M_PI)*sigma);
        f += w[0]*4.0/(1.0 + pow(x(0), 2.0));
    }

    std::cout << f << std::endl;

	return 0;
}