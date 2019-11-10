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
			
	//----------Save initial value----------
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
			LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
			std::vector<double> Q = std::vector<double>(KDEGREE, 0.0);		//Interna load vector
			std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
			for (auto element : elements) {
				Matrix<double> Ke;
				Vector<double> Qe;
				TotalLagrangeSolid<double, ShapeFunction8Cubic, Gauss8Cubic>(Ke, Qe, nodes, u, element, 1000.0, 0.3);
				Assembling(K, Q, Ke, Qe, element, field);
			}

			//----------Set Neumann Boundary Condition----------
			SetNeumann(F, isqfixed, qfixed);

			//----------Get residual load vector----------
			std::vector<double> R = std::vector<double>(KDEGREE, 0.0);		//Residual load vector
			for(int i = 0; i < KDEGREE; i++){
				R[i] = (finc / (double)fincmax)*F[i] - Q[i];
			}
						
			//----------Set Dirichlet Boundary Condition----------
			SetDirichlet(K, R, isufixed, ufixed, 1.0e10);

			//----------Check convergence----------
			double normR = 0.0, normF = 0.0;
			for(int i = 0; i < KDEGREE; i++){
				normR += pow(R[i], 2.0);
				normF += pow((finc / (double)fincmax)*F[i], 2.0);
			}
			normR = sqrt(normR);
			normF = sqrt(normF);
			if (normR / normF < 1.0e-10) {
				std::cout << "k = " << k << "\tNorm = " << normR / normF << std::endl;
				break;
			}
			
			//----------Solve System Equation----------
			CSR<double> Kmod = CSR<double>(K);
			std::vector<double> result = ScalingCG(Kmod, R, 100000, 1.0e-10);

			//----------Update displacement u----------
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