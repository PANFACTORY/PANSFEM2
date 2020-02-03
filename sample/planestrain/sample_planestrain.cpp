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
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/planestrain/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");
	
	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Add Neumann Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, model_path + "Neumann.csv");

	//----------Culculate Ke Fe and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for (auto element : elements) {
		Matrix<double> Ke;
		PlaneStrain<double, ShapeFunction3Triangle, Gauss1Triangle>(Ke, nodes, element, 210000.0, 0.3, 1.0);
		Assembling(K, Ke, element, field);
	}

	//----------Set Neumann Boundary Condition----------
	SetNeumann(F, isqfixed, qfixed);

	//----------Set Dirichlet Boundary Condition----------
	SetDirichlet(K, F, isufixed, ufixed, 1.0e10);

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
	std::vector<int> et = std::vector<int>(elements.size(), 5);
	AddElementTypes(et, fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	Vector<double> Fe;
	std::vector<int> Se = { 3, 4, 5 };
	PlaneSurfaceForce<double, ShapeFunction3Line, Gauss2Line>(Fe, nodes, Se, 1.0, 1.0, 1.0);
	std::cout << Fe << std::endl;


	return 0;
}