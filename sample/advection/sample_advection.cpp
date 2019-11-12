#include <iostream>
#include <vector>
#include <cmath>


#include "LinearAlgebra/Models/Vector.h"
#include "LinearAlgebra/Models/Matrix.h"
#include "LinearAlgebra/Models/LILCSR.h"
#include "PrePost/Import/ImportFromCSV.h"
#include "FEM/Controller/Assembling.h"
#include "FEM/Equation/Advection.h"
#include "FEM/Controller/BoundaryCondition.h"
#include "LinearAlgebra/Solvers/CG.h"
#include "PrePost/Export/ExportToVTK.h"
#include "FEM/Controller/ShapeFunction.h"
#include "FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/advection/";
	
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

	//----------Initialize T----------
	std::vector<double> T = std::vector<double>(nodes.size(), 0.0);
	T[12] = 100.0;

	//----------Define time step and theta----------
	double dt = 0.1;
	double theta = 0.5;

	//----------Time step loop----------
	for(int t = 0; t < 1000; t++){
		std::cout << "t = " << t << std::endl;

		//----------Culculate Ke Fe and Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
		for (auto element : elements) {
			std::vector<double> Ttmp;
			for(auto i : element){
				Ttmp.push_back(T[i]);
			}
			Vector<double> Te = Vector<double>(Ttmp);

			Matrix<double> Ke;
			Advection<double, ShapeFunction3Triangle, Gauss1Triangle>(Ke, nodes, element, 1.0, 1.0, 1.0);
			Matrix<double> Ce;
			Mass<double, ShapeFunction3Triangle, Gauss1Triangle>(Ce, nodes, element, 1.0);
			Matrix<double> Ae = Ce / dt + Ke * theta;
			Vector<double> be = (Ce / dt - Ke * (1.0 - theta)) * Te; 
			Assembling(K, F, Ae, be, element, field);
		}

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(K, F, isufixed, ufixed, 1.0e4);

		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(K);
		CSR<double> M = ILU0(Kmod);
		std::vector<double> result = ILU0BiCGSTAB(Kmod, M, F, 10000, 1.0e-10);

		//----------Update T----------
		T.clear();
		FieldResultToNodeValue(result, T, field);
				
		//----------Save initial value----------
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(elements.size(), 5);
		AddElementTypes(et, fout);
		AddPointScalers(T, "T", fout);
		fout.close();
	}

	return 0;
}