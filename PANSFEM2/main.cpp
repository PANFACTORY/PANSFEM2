#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/ImportExport/Import.h"
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

	//----------�S�́\�ߓ_�֌W----------
	

	//----------�A�Z���u�����O----------
	LILCSR<double> K = LILCSR<double>(12, 12);
	std::vector<double> F = std::vector<double>(12, 0.0);
	for (auto element : elements) {
		PlaneStrainTri(K, nodes, element, 210000.0, 0.3, 1.0);
	}

	//----------Dirichlet���E�����̓K�p----------
	std::vector<int> isufixed = { 0, 1, 10 };
	std::vector<double> ufixed = { 0.0, 0.0, 0.0 };
	SetDirichlet(K, F, isufixed, ufixed, 1.0e20);

	//----------Neumann���E�����̓K�p----------
	std::vector<int> isqfixed = { 7 };
	std::vector<double> qfixed = { -100.0 };
	SetNeumann(F, isqfixed, qfixed);

	//----------�S�́\�ߓ_������������----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> u = CG(Kmod, F, 100, 1.0e-10);

	//----------���̏o��----------
	for (auto ui : u) {
		std::cout << ui << std::endl;
	}
	
	return 0;
}