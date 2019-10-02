#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/Point.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "FEM/Controller/PlaneStrain.h"
#include "FEM/Controller/Dirichlet.h"
#include "LinearAlgebra/Solvers/CG.h"


using namespace PANSFEM2;


int main() {
	//----------�ߓ_�ǉ�----------
	std::vector<Point<double> > nodes;
	nodes.push_back(Point<double>(0.0, 0.0));
	nodes.push_back(Point<double>(1.0, 0.0));
	nodes.push_back(Point<double>(2.0, 0.0));
	nodes.push_back(Point<double>(2.0, 1.0));
	nodes.push_back(Point<double>(1.0, 1.0));
	nodes.push_back(Point<double>(0.0, 1.0));
	
	//----------�v�f�\�ߓ_�֌W----------
	std::vector<std::vector<int> > nodestoelements;
	nodestoelements.push_back({ 0, 1, 4 });
	nodestoelements.push_back({ 1, 2, 3 });
	nodestoelements.push_back({ 3, 4, 1 });
	nodestoelements.push_back({ 4, 5, 0 });

	//----------�A�Z���u�����O----------
	LILCSR<double> K = LILCSR<double>(12, 12);
	std::vector<double> F = std::vector<double>(12, 0.0);
	for (auto nodestoelement : nodestoelements) {
		PlaneStrainTri(K, F, nodes, nodestoelement, 210000.0, 0.3, 1.0);
	}

	//----------Dirichlet���E�����̓K�p----------
	std::vector<int> isufixed = { 0, 1, 10 };
	std::vector<double> ufixed = { 0.0, 0.0, 0.0 };
	SetDirichlet(K, F, isufixed, ufixed, 1.0e20);

	//----------Neumann���E�����̓K�p----------
	F[7] = -100.0;

	//----------�S�́\�ߓ_������������----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> u = CG(Kmod, F, 100, 1.0e-10);

	//----------���̏o��----------
	for (auto ui : u) {
		std::cout << ui << std::endl;
	}

	return 0;
}