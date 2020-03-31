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
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/NewtonCotesIntegration.h"
#include "../../src/PrePost/Import/ImportFromVTK.h"


using namespace PANSFEM2;


Vector<double> P(Vector<double> _x){
	double sigma = 1.0;
	double mu = 300.0;
	double dmu = 0.0;
    return { 0.0, -exp(-0.5*pow((_x(0) - (mu + dmu))/sigma, 2.0))/(sqrt(2.0*M_PI)*sigma) };
}


int main() {
	//----------Model Path----------
	std::string model_path = "sample/concentrationgauss/";
	
	//----------Add Nodes----------
    ImportModelFromVTK<double> modelimport(model_path + "result223.vtk", 2);
	std::vector<Vector<double> > nodes = modelimport.ImportPOINTS();
	std::vector<std::vector<int> > elements = modelimport.ImportCELLS();

    //----------Add Edge----------
	std::vector<std::vector<int> > edges;
	ImportElementsFromCSV(edges, model_path + "Edge.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 823.0;
    double E2 = 210000.0;
	double Poisson = 0.3;
    double rho0 = 0.0;
    double rho1 = 0.0323;
    double rho2 = 1.0;
	double p = 3.0;
    double q = 3.0;

    std::vector<double> r0 = modelimport.ImportCELLSCALARS("r0");
    std::vector<double> r1 = modelimport.ImportCELLSCALARS("r1");
        
    //*************************************************
    //  Get robust compliance value and sensitivities
    //*************************************************
    //----------Assembling----------
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
    for (int i = 0; i < elements.size(); i++) {
        double E = E0*(1.0 - pow(r0[i], p)) + (E1*(1.0 - pow(r1[i], q)) + E2*pow(r1[i], q))*pow(r0[i], p);
        Matrix<double> Ke;
        PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
        Assembling(K, Ke, elements[i], field);
    }

    //----------Culculate Surface force and Assembling----------
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for(auto edge : edges){
		Vector<double> Fe;
		PlaneSurfaceForce<double, ShapeFunction2Line, NewtonCotes3Line>(Fe, nodes, edge, P, 1.0);
		Assembling(F, Fe, edge, field);
	}

    //----------Set Dirichlet Boundary Condition(Only fix displacement 0.0)----------
    SetDirichlet(K, isufixed, ufixed, 1.0e10);
    CSR<double> Kmod = CSR<double>(K);	

    //----------Get function value and sensitivities----------
    std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);
    double f = std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
	
    std::cout << f << std::endl;
	
	return 0;
}