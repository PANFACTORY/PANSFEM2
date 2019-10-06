#include <iostream>
#include <vector>


#include "../../../LinearAlgebra/Models/Vector.h"
#include "../../../LinearAlgebra/Models/LILCSR.h"
#include "../../../PrePost/Import/ImportFromCSV.h"
#include "../../../FEM/Controller/Assembling.h"
#include "../../../FEM/Equation/HeatTransfer.h"
#include "../../../FEM/Controller/BoundaryCondition.h"
#include "../../../LinearAlgebra/Solvers/CG.h"
#include "../../../PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


void HeatTransferSample() {
	//----------Model Path----------
	std::string model_path = "Samples/HeatTransfer/Static/";

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

	//----------Add Initial Condition----------
	std::vector<double> T = std::vector<double>(nodes.size(), 0.0);

	//----------Time Step----------
	double dt = 0.1;
	double theta = 0.5;
	for (int k = 0; k < 100; k++) {
		//----------Culculate Ke Ce and Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		LILCSR<double> C = LILCSR<double>(KDEGREE, KDEGREE);
		for (auto element : elements) {
			std::vector<std::vector<double> > Ke = HeatTransferTri(nodes, element, 1.0, 1.0);
			Assembling(K, Ke, element, field);
			std::vector<std::vector<double> > Ce = HeatCapacityTri(nodes, element, 1.0, 1.0, 1.0);
			Assembling(C, Ce, element, field);
		}

		//----------Make Equation----------
		LILCSR<double> KC = (1.0 / dt)*C + theta * K;
		std::vector<double> F = ((1.0 / dt)*C - (1.0 - theta) * K) * T;

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(KC, F, isufixed, ufixed, 1.0e3);

		//----------Set Neumann Boundary Condition----------
		SetNeumann(F, isqfixed, qfixed);

		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(KC);
		std::vector<double> result = CG(Kmod, F, 100, 1.0e-10);

		//----------Post Process----------
		T.clear();
		FieldResultToNodeValue(result, T, field);

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(32, 5);
		AddElementTypes(et, fout);
		AddPointScalers(T, "T", fout);
		fout.close();
	}
}