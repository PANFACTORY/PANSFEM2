#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/ImportExport/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/PlaneStrain.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"


using namespace PANSFEM2;


int main() {
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, "Data/Solid/TriangleBeam/Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, "Data/Solid/TriangleBeam/Element.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), "Data/Solid/TriangleBeam/Field.csv");

	//----------Culculate Ke and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
	for (auto element : elements) {
		std::vector<std::vector<double> > Ke = PlaneStrainTri(nodes, element, 210000.0, 0.3, 1.0);
		Assembling(K, Ke, element, field);
	}
	std::cout << K << std::endl;

	//----------Set Dirichlet Boundary Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, "Data/Solid/TriangleBeam/Dirichlet.csv");
	SetDirichlet(K, F, isufixed, ufixed, 1.0e20);

	//----------Set Neumann Boundary Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, "Data/Solid/TriangleBeam/Neumann.csv");
	SetNeumann(F, isqfixed, qfixed);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = CG(Kmod, F, 100, 1.0e-10);

	//----------Post Process----------
	std::vector<Vector<double> > u;
	FieldResultToNodeValue(result, u, field);
	for (auto ui : u) {
		std::cout << ui;
	}
	
	return 0;
}