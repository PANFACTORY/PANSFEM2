#include <iostream>
#include <vector>
#include <cmath>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/Matrix.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/Import/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/Solid.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"
#include "PrePost/Export/ExportToVTK.h"
#include "FEM/Controller/ShapeFunction.h"
#include "FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/Optimize/";
	
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

	//----------Initialize u----------
	std::vector<Vector<double> > u;
	
	//----------Culculate Ke Fe and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for (auto element : elements) {
		Matrix<double> Ke;
		LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic>(Ke, nodes, element, 210000.0, 0.3);
		Assembling(K, Ke, element, field);
	}

	//----------Set Dirichlet Boundary Condition----------
	SetDirichlet(K, F, isufixed, ufixed, 1.0e10);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = ScalingCG(Kmod, F, KDEGREE, 1.0e-10);

	//----------Update T----------
	u.clear();
	FieldResultToNodeValue(result, u, field);
			
	//----------Save initial value----------
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	std::vector<int> et = std::vector<int>(elements.size(), 25);
	AddElementTypes(et, fout);
	AddPointVectors(u, "u", fout);
	fout.close();
	
	return 0;
}