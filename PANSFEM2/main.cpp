#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/Import/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/TotalLagrange.h"
#include "FEM/Equation/Solid.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"
#include "PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "Samples/TotalLagrange/";

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

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Add Neumann Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, model_path + "Neumann.csv");

	//----------Add Initial Condition----------
	std::vector<std::vector<double> > u = std::vector<std::vector<double> >(nodes.size(), std::vector<double>(3, 0.0));
	
	//----------Save file----------
	std::ofstream fout(model_path + "result" + std::to_string(0) + ".vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	std::vector<int> et = std::vector<int>(1250, 12);
	AddElementTypes(et, fout);
	AddPointVectors(u, "u", fout);
	fout.close();

	//----------Inclement Load----------
	for (int finc = 1, fincmax = 100; finc <= fincmax; finc++) {
		std::cout << "finc = " << finc << "\t";

		//----------Newton-Raphson Step----------
		for (int k = 0; k < 10; k++) {
			//----------Culculate Ke Fe and Assembling----------
			LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
			std::vector<double> Q = std::vector<double>(KDEGREE, 0.0);
			std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
			for (auto element : elements) {
				std::vector<std::vector<double> > Ke;
				std::vector<double> Qe;
				TotalLagrange(Ke, Qe, nodes, u, element, 1000.0, 0.3);
				Assembling(K, Q, Ke, Qe, element, field);
			}

			//----------Set Neumann Boundary Condition----------
			SetNeumann(F, isqfixed, qfixed);
			std::vector<double> R = (finc / (double)fincmax)*F - Q;
			
			//----------Set Dirichlet Boundary Condition----------
			SetDirichlet(K, R, isufixed, ufixed, 1.0e10);
			if (Norm(R) / Norm((finc / (double)fincmax)*F) < 1.0e-10) {
				std::cout << "k = " << k << "\tNorm = " << Norm(R) / Norm((finc / (double)fincmax)*F) << std::endl;
				break;
			}
			
			//----------Solve System Equation----------
			CSR<double> Kmod = CSR<double>(K);
			//CSR<double> M = ILU0(Kmod);
			std::vector<double> result = ScalingCG(Kmod, R, 10000, 1.0e-10);

			//----------Post Process----------
			std::vector<std::vector<double> > du;
			FieldResultToNodeValue(result, du, field);
			u += du;
		}

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(finc) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(et, fout);
		AddPointVectors(u, "u", fout);
		fout.close();
	}
	
	return 0;
}