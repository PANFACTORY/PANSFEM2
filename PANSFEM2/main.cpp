#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/ImportExport/Import.h"
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

	//----------ëSëÃÅ\êﬂì_ä÷åW----------
	std::vector<std::pair<int, int> > systemindex;
	systemindex.push_back(std::make_pair(0, 2));
	systemindex.push_back(std::make_pair(2, 2));
	systemindex.push_back(std::make_pair(4, 2));
	systemindex.push_back(std::make_pair(6, 2));
	systemindex.push_back(std::make_pair(8, 2));
	systemindex.push_back(std::make_pair(10, 2));

	//----------Culculate Ke and Assembling----------
	LILCSR<double> K = LILCSR<double>(12, 12);
	std::vector<double> F = std::vector<double>(12, 0.0);
	for (auto element : elements) {
		std::vector<std::vector<double> > Ke;
		PlaneStrainTri(Ke, nodes, element, 210000.0, 0.3, 1.0);
		Assembling(K, Ke, element, systemindex);
	}

	std::cout << K << std::endl;

	//----------Set Dirichlet Boundary Condition----------
	std::vector<int> isufixed = { 0, 1, 10 };
	std::vector<double> ufixed = { 0.0, 0.0, 0.0 };
	SetDirichlet(K, F, isufixed, ufixed, 1.0e20);

	//----------Set Neumann Boundary Condition----------
	std::vector<int> isqfixed = { 7 };
	std::vector<double> qfixed = { -100.0 };
	SetNeumann(F, isqfixed, qfixed);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> u = CG(Kmod, F, 100, 1.0e-10);

	//----------Post Process----------
	for (auto ui : u) {
		std::cout << ui << std::endl;
	}
	
	return 0;
}