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

	//----------Initialize design variables----------
	std::vector<double> xk = std::vector<double>(elements.size(), 0.5);
	std::vector<double> xkm1 = std::vector<double>(elements.size());
	std::vector<double> xkm2 = std::vector<double>(elements.size());

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	std::vector<double> L = std::vector<double>(elements.size());
	std::vector<double> U = std::vector<double>(elements.size());
	double s = 0.7;
	
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
			double E = E1 * pow(xk[i], p) + E0 * (1.0 - pow(xk[i], p));
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
		AddElementScalers(xk, "s", fout);
		fout.close();


		//**************************************************
		//	Update design variables with MMA
		//**************************************************
		double compliance = 0.0;																				//Function value of compliance
		std::vector<double> p0 = std::vector<double>(elements.size(), 0.0);										//Positive sensitivities of compliance
		std::vector<double> q0 = std::vector<double>(elements.size(), 0.0);										//Negative sensitivities of compliance
		double weight = 0.0;																					//Function values of weight
		std::vector<Vector<double> > ps = std::vector<Vector<double> >(elements.size(), Vector<double>(1));		//Positive sensitivities of weight
		std::vector<Vector<double> > qs = std::vector<Vector<double> >(elements.size(), Vector<double>(1));		//Negative sensitivities of weight

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
			compliance += (E1 * pow(xk[i], p) + E0 * (1.0 - pow(xk[i], p))) * ueKeue;
			double dcompliance = p * (E1 * pow(xk[i], p - 1.0) - E0 * pow(xk[i], p - 1.0)) * ueKeue;
			if (dcompliance > 0.0) {
				p0[i] = dcompliance;
			} else {
				q0[i] = dcompliance;
			}

			//.....Constraint functions.....
			weight += xk[i] - weightlimit;
			double dweight = 1.0;
			if (dweight > 0.0) {
				ps[i](0) = dweight;
			} else {
				qs[i](0) = dweight;
			}
		}

		//----------Check convergence----------
		if(fabs((compliance - compliancebefore) / (compliance + compliancebefore)) < complianceeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}

		std::cout << "Compliance:\t" << compliance << "\t";
		std::cout << "Weight:\t" << weight << "\t";

		//----------Set parameter U and L----------
		if(k < 2){
			for(int i = 0; i < elements.size(); i++){
				L[i] = xk[i] - (1.0 - 0.0);
				U[i] = xk[i] + (1.0 - 0.0);
			}
		} else {
			for(int i = 0; i < elements.size(); i++){
				if((xk[i] - xkm1[i])*(xkm1[i] - xkm2[i]) < 0.0){
					L[i] = xk[i] - s*(xkm1[i] - L[i]);
					U[i] = xk[i] + s*(U[i] - xkm1[i]);
				} else {
					L[i] = xk[i] - (xkm1[i] - L[i])/s;
					U[i] = xk[i] + (U[i] - xkm1[i])/s;
				}
			}
		}

		//----------Set move limit----------
		std::vector<double> xmin = std::vector<double>(elements.size());
		std::vector<double> xmax = std::vector<double>(elements.size());
		for(int i = 0; i < elements.size(); i++){
			xmin[i] = 0.9*L[i] + 0.1*xk[i];
			xmax[i] = 0.9*U[i] + 0.1*xk[i];
		}

		//----------Loop for solving subproblem----------
		std::vector<double> xkp1 = std::vector<double>(elements.size());
		Vector<double> y = Vector<double>(1);
		for(int t = 0; t < 100; t++){
			//.....Get x(y).....
			for(int i = 0; i < elements.size(); i++){
				double dlmin = (p0[i] + y*ps[i]) / pow(U[i] - xmin[i], 2.0) - (q0[i] + y*qs[i]) / pow(xmin[i] - L[i], 2.0);
				double dlmax = (p0[i] + y*ps[i]) / pow(U[i] - xmax[i], 2.0) - (q0[i] + y*qs[i]) / pow(xmax[i] - L[i], 2.0);

				if (dlmax >= 0.0) {
					xkp1[i] = xmin[i];
				} else if (dlmax <= 0.0) {
					xkp1[i] = xmax[i];
				} else {
					xkp1[i] = (sqrt(p0[i] + y*ps[i])*L[i] + sqrt(q0[i] + y*qs[i])*U[i]) / (sqrt(p0[i] + y*ps[i]) + sqrt(q0[i] + y*qs[i]));
				}
			}

			//.....Get y.....

		}

		//----------Update xk----------
		

		
	}
	
	return 0;
}