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
	std::string model_path = "sample/Optimize/";
	
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
	
	//----------Optimize loop----------
	for(int i = 0; i < 100; i++){
		//**************************************************
		//	Excute direct analysis
		//**************************************************

		//----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
		for (auto element : elements) {
			Matrix<double> Ke;
			LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic >(Ke, nodes, element, 210000.0, 0.3);
			Assembling(K, Ke, element, field);
		}
		
		//----------Set Neumann Boundary Condition----------
		SetNeumann(F, isqfixed, qfixed);
		
		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(K, F, isufixed, ufixed, 1.0e10);
		
		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
		
		//----------Post Process----------
		std::vector<Vector<double> > u;
		FieldResultToNodeValue(result, u, field);

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(i) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(elements.size(), 25);
		AddElementTypes(et, fout);
		AddPointVectors(u, "u", fout);
		fout.close();


		//**************************************************
		//	Get sensitivity and update design variables
		//**************************************************

		
	}
	
	return 0;
}