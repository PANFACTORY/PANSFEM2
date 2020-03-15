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


using namespace PANSFEM2;


Vector<double> f(Vector<double> _x){
	double sigma = 10.0;
	double mu = 30.0;
	double xmax = 100.0;
	double xmin = 0.0;
	double delta = 0.5*(erf((xmax - mu)/(sqrt(2.0)*sigma))- erf((xmin - mu)/(sqrt(2.0)*sigma)));
    return { 0.0, -exp(-0.5*pow((_x(0) - mu)/sigma, 2.0))/(sqrt(2.0*M_PI)*sigma)/delta };
}


int main() {
	//----------Model Path----------
	std::string model_path = "sample/concentrationgauss/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

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

	//----------Culculate Ke, body force and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for (auto element : elements) {
		Matrix<double> Ke;
		PlaneStrain<double, ShapeFunction4Square, Gauss4Square>(Ke, nodes, element, 210000.0, 0.3, 1.0);
		Assembling(K, Ke, element, field);
	}

	//----------Culculate Surface force and Assembling----------
	for(auto edge : edges){
		Vector<double> Fe;
		PlaneSurfaceForce<double, ShapeFunction2Line, NewtonCotes3Line>(Fe, nodes, edge, f, 1.0);
		Assembling(F, Fe, edge, field);
	}

	//----------Set Dirichlet Boundary Condition----------
	SetDirichlet(K, F, isufixed, ufixed, 1.0e9);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);

	//----------Update displacement u----------
	std::vector<Vector<double> > u = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
	FieldResultToNodeValue(result, u, field);
			
	//----------Save initial value----------
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 9), fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	for(int i = 10; i < nodes.size(); i+=11){
		std::cout << F[2*i] << "\t" << F[2*i + 1] << std::endl;
	}

	return 0;
}