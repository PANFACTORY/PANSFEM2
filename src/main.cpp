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

	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	std::vector<double> L = std::vector<double>(elements.size());
	std::vector<double> U = std::vector<double>(elements.size());
	
	double weightlimit = 0.5;
	double compliancebefore = 0.0;
	double complianceeps = 1.0e-5;
	
	//----------Optimize loop----------
	for(int k = 0; k < 50; k++){
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
			LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic >(Ke, nodes, elements[i], E, Poisson);
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
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(elements.size(), 25);
		AddElementTypes(et, fout);
		AddPointVectors(u, "u", fout);
		AddElementScalers(s, "s", fout);
		fout.close();


		//**************************************************
		//	Get sensitivity and update design variables
		//**************************************************
		double compliance = 0.0;													//Function value of compliance
		std::vector<double> dcompliances = std::vector<double>(elements.size());	//Sensitivities of compliance
		double weight = 0.0;														//Function values of weight
		std::vector<double> dweights = std::vector<double>(elements.size());		//Sensitivities of weight

		//----------Get function values and sensitivities----------
		for (int i = 0; i < elements.size(); i++) {
			Vector<double> ue = Vector<double>();
			for(int j = 0; j < elements[i].size(); j++){
				ue = ue.Vstack(u[elements[i][j]]);
			}

			Matrix<double> Ke;
			LinearIsotropicElasticSolid<double, ShapeFunction20Cubic, Gauss27Cubic >(Ke, nodes, elements[i], 1.0, Poisson);
			double ueKeue = (ue.Transpose()*Ke*ue)(0);

			compliance += (E1 * pow(s[i], p) + E0 * (1.0 - pow(s[i], p))) * ueKeue;
			dcompliances[i] += p * (E1 * pow(s[i], p - 1.0) - E0 * pow(s[i], p - 1.0)) * ueKeue;
			weight += s[i] - weightlimit;
			dweights[i] = 1.0;
		}

		//----------Check convergence----------
		if(fabs((compliance - compliancebefore) / (compliance + compliancebefore)) < complianceeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}

		std::cout << "Compliance:\t" << compliance << "\t";
		std::cout << "Weight:\t" << weight << "\t";

		//----------Get updated design variables with MMA----------
		//----------Set parameter U and L----------
		if(k < 2){

		} else {

		}

		//----------Set move limit----------
		std::vector<double> smin = std::vector<double>(elements.size());
		std::vector<double> smax = std::vector<double>(elements.size());
		for(int i = 0; i < elements.size(); i++){
			
		}



		double lambda0 = lambdamin, lambda1 = lambdamax, lambda;
		std::vector<double> snext = std::vector<double>(elements.size());
		while((lambda1 - lambda0) / (lambda1 + lambda0) > lambdaeps){
			lambda = 0.5 * (lambda1 + lambda0);

			for (int i = 0; i < elements.size(); i++) {
				snext[i] = pow(dcompliances[i] / (dweights[i] * lambda), iota) * s[i];
				if(snext[i] < std::max(0.0, (1.0 - movelimit)*s[i])) {
					snext[i] = std::max(0.0, (1.0 - movelimit)*s[i]);
				} else if(snext[i] > std::min(1.0, (1.0 + movelimit)*s[i])) {
					snext[i] = std::min(1.0, (1.0 + movelimit)*s[i]);
				}
			}

			double weightnext = 0.0;
			for (int i = 0; i < elements.size(); i++) {
				weightnext += snext[i] - weightlimit;
			}

			if (weightnext > 0.0) {
				lambda0 = lambda;
			}
			else {
				lambda1 = lambda;
			}
		}

		std::cout << "Lagrange value:\t" << lambda << std::endl;

		//----------Update design variables and compliance----------
		for (int i = 0; i < elements.size(); i++) {
			s[i] = snext[i];
		}
		compliancebefore = compliance;
	}
	
	return 0;
}