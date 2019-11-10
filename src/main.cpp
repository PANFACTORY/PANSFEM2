#include <iostream>
#include <vector>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/Matrix.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/Import/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/Solid.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"
#include "PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/solid/";
	
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

	//----------Set initial displacement u----------
	std::vector<Vector<double> > u = std::vector<Vector<double> >(nodes.size(), Vector<double>(3)); 
			
	//----------Save file----------
	std::ofstream fout(model_path + "result" + std::to_string(0) + ".vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	std::vector<int> et = std::vector<int>(elements.size(), 12);
	AddElementTypes(et, fout);
	AddPointVectors(u, "u", fout);
	fout.close();

	//----------Inclement Load----------
	for (int finc = 1, fincmax = 100; finc <= fincmax; finc++) {
		std::cout << "finc = " << finc << "\t";

		//----------Newton-Raphson Step----------
		for (int k = 0; k < 100; k++) {
			//----------Culculate Ke Fe and Assembling----------
			LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
			std::vector<double> Q = std::vector<double>(KDEGREE, 0.0);
			std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
			for (auto element : elements) {
				Matrix<double> Ke;
				Vector<double> Qe;
				TotalLagrangeSolid(Ke, Qe, nodes, u, element, 1000.0, 0.3);
				Assembling(K, Q, Ke, Qe, element, field);
			}

			//----------Set Neumann Boundary Condition----------
			SetNeumann(F, isqfixed, qfixed);
			std::vector<double> R = (finc / (double)fincmax)*F - Q;
			
			//----------Set Dirichlet Boundary Condition----------
			SetDirichlet(K, R, isufixed, ufixed, 1.0e10);

			//----------Check convergence----------
			if (Norm(R) / Norm((finc / (double)fincmax)*F) < 1.0e-10) {
				std::cout << "k = " << k << "\tNorm = " << Norm(R) / Norm((finc / (double)fincmax)*F) << std::endl;
				break;
			}
			std::cout << "\tNorm = " << Norm(R) / Norm((finc / (double)fincmax)*F) << std::endl;
			
			//----------Solve System Equation----------
			CSR<double> Kmod = CSR<double>(K);
			std::vector<double> result = ScalingCG(Kmod, R, 100000, 1.0e-10);

			//----------Post Process----------
			std::vector<Vector<double> > du;
			FieldResultToNodeValue(result, du, field);
			for(int i = 0; i < nodes.size(); i++){
				u[i] += du[i];
			}
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