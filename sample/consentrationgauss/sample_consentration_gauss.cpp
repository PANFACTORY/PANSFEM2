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
    double f = 0.0;
    for(auto element : elements){
        Matrix<double> xe = Matrix<double>(0, 1);
		for(auto i : element){
			xe= xe.Vstack(x[i].Transpose());
		}

        for(int g = 0; g < NewtonCotes3Line<double>::N; g++){
            Vector<double> N = ShapeFunction2Line<double>::N(NewtonCotes3Line<double>::Points[g]);
            Matrix<double> dNdr = ShapeFunction2Line<double>::dNdr(NewtonCotes3Line<double>::Points[g]);

			//----------Get difference of shape function----------
			Vector<double> X = xe.Transpose()*N;
            Matrix<double> dXdr = dNdr*xe;
			double dl = sqrt((dXdr*dXdr.Transpose())(0, 0));
            double df = 4.0/(1.0 + pow(X(0), 2.0));
            f += NewtonCotes3Line<double>::Weights[g][0]*df*dl;
        }
    }

    std::cout << f << std::endl;
	return 0;
}