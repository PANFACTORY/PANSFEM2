/*#pragma once
#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/Import/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/PlaneStrain.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"
#include "PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "Samples/Solid/TriangleBeam/";

	//----------Add Nodes----------
	std::vector<std::vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");

	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

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
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");
	SetDirichlet(K, F, isufixed, ufixed, 1.0e20);

	//----------Set Neumann Boundary Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, model_path + "Neumann.csv");
	SetNeumann(F, isqfixed, qfixed);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = CG(Kmod, F, 100, 1.0e-10);

	//----------Post Process----------
	std::vector<std::vector<double> > u;
	FieldResultToNodeValue(result, u, field);

	//----------Save file----------
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes({ 5, 5, 5, 5 }, fout);
	AddPointVectors(u, "u", fout);
	fout.close();

	return 0;
}*/