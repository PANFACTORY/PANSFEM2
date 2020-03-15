#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Beam.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/beam/";
	
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

    //----------Generate director vector----------
    std::vector<Vector<double> > v1 = std::vector<Vector<double> >(nodes.size(), { 0.0, 1.0, 0.0 });
	std::vector<Vector<double> > v2 = std::vector<Vector<double> >(nodes.size(), { 0.0, 0.0, 1.0 });

	//----------Culculate Ke, body force and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for (auto element : elements) {
		Matrix<double> Ke;
		LinearIsotropicElasticBeam<double, ShapeFunction2Line, Gauss1Line, Gauss4Square>(Ke, nodes, element, v1, v2, 210000.0, 0.3, 1.0, 1.0);
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
	std::vector<Vector<double> > ur = std::vector<Vector<double> >(nodes.size(), Vector<double>(6));
	FieldResultToNodeValue(result, ur, field);
    std::vector<Vector<double> > u = std::vector<Vector<double> >(nodes.size());
    for(int i = 0; i < nodes.size(); i++){
        u[i] = ur[i].Segment(0, 3);
    }
			
	//----------Save initial value----------
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	std::vector<int> et = std::vector<int>(elements.size(), 3);
	AddElementTypes(et, fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	return 0;
}