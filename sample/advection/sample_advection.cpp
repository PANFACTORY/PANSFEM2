#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Advection.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


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
	T[5000] = 1.0;	T[5101] = 1.0;	T[5202] = 1.0;
	T[4999] = 1.0;	T[5100] = 2.0;	T[5201] = 1.0;
	T[4998] = 1.0;	T[5099] = 1.0;	T[5200] = 1.0;
	
	//----------Define time step and theta----------
	double dt = 0.001;
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
			Advection<double, ShapeFunction4Square, Gauss4Square>(Ke, nodes, element, 0.1, 0.1, 1.0);
			Matrix<double> Ce;
			Mass<double, ShapeFunction4Square, Gauss4Square>(Ce, nodes, element, 1.0);
			Matrix<double> Ae = Ce/dt + theta*Ke;
			Vector<double> be = (Ce/dt - (1.0 - theta)*Ke)*Te; 
			Assembling(K, F, Ae, be, element, field);
		}

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(K, F, isufixed, ufixed, 1.0e5);

		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = ScalingBiCGSTAB(Kmod, F, 100000, 1.0e-10);

		//----------Update T----------
		FieldResultToNodeValue(result, T, field);
				
		//----------Save initial value----------
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointScalers(T, "T", fout, true);
		fout.close();
	}

	return 0;
}