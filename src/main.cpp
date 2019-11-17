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
#include "Optimize/Solver/MMA.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/optimize/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > X;
	ImportNodesFromCSV(X, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, X.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Add Neumann Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, model_path + "Neumann.csv");

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;
	double weightlimit = 0.5;

	//----------Initialize optimization solver----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);
	MMA<double> optimizer = MMA<double>(elements.size(), 1);
	
	//----------Optimize loop----------
	for(int k = 0; k < 5; k++){
		std::cout << "k = " << k << "\t";

		//**************************************************
		//	Excute direct analysis
		//**************************************************

		//----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1 * pow(s[i], p) + E0 * (1.0 - pow(s[i], p));
			Matrix<double> Ke;
			LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic >(Ke, X, elements[i], E, Poisson);
			Assembling(K, Ke, elements[i], field);
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
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(X, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(elements.size(), 25);
		AddElementTypes(et, fout);
		AddPointVectors(u, "u", fout);
		AddElementScalers(s, "s", fout);
		fout.close();


		//**************************************************
		//	Update design variables
		//**************************************************
		double compliance = 0.0;																						//Function value of compliance
		std::vector<double> dcompliance = std::vector<double>(elements.size());											//Sensitivities of compliance
		Vector<double> constraints = Vector<double>(1);																	//Function values of weight
		std::vector<Vector<double> > dconstraints = std::vector<Vector<double> >(elements.size(), Vector<double>(1));	//Sensitivities of weight

		//----------Get function values and sensitivities----------
		for (int i = 0; i < elements.size(); i++) {
			//.....Objective function.....
			Vector<double> ue = Vector<double>();
			for(int j = 0; j < elements[i].size(); j++){
				ue = ue.Vstack(u[elements[i][j]]);
			}
			Matrix<double> Ke;
			LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic >(Ke, X, elements[i], 1.0, Poisson);
			double ueKeue = (ue.Transpose()*Ke*ue)(0);
			compliance += (E1 * pow(s[i], p) + E0 * (1.0 - pow(s[i], p))) * ueKeue;
			dcompliance[i] = -p * (E1 * pow(s[i], p - 1.0) - E0 * pow(s[i], p - 1.0)) * ueKeue;
			
			//.....Constraint functions.....
			constraints(0) += s[i] - weightlimit;
			dconstraints[i](0) = 1.0;
		}

		std::cout << "Compliance:\t" << compliance << "\t";
		std::cout << "Weight:\t" << constraints(0) << "\t";

		//----------Check convergence----------
		if(optimizer.IsConvergence(compliance)){
			std::cout << "--------------------Optimized--------------------" << std::endl;
			break;
		}

		//----------Update s----------
		optimizer.UpdateVariables(s, compliance, dcompliance, constraints, dconstraints);	
	}

	return 0;
}