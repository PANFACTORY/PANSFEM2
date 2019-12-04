#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/Lanczos.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/optimize_bibration/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, x.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.0001;
	double E1 = 100000.0;
	double Poisson = 0.3;
	double p = 3.0;

    double rho0 = 1.0e-12;
    double rho1 = 1.0e-8;

    int m = 5;

	int n = 1;

	double iota = 0.75;
	double lambdamin = 1.0e-20;
	double lambdamax = 1.0e20;
	double lambdaeps = 1.0e-15;
	double movelimit = 0.15;

	double weightlimit = 0.32;
	double objectivebefore = 0.0;
	double objectiveeps = 1.0e-5;
	
	//----------Optimize loop----------
	for(int k = 0; k < 200; k++){
		std::cout << "k = " << k << "\t";

		//**************************************************
		//	Excute direct analysis
		//**************************************************

		//----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		LILCSR<double> M = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1 * pow(s[i], p) + E0 * (1.0 - pow(s[i], p));
            double rho = rho1*s[i] + rho0*(1.0 - s[i]); 
			if(i == 1389 || i == 1390 || i == 1409 || i == 1410){
				rho += 1.25e-6;
			}
			Matrix<double> Ke, Me;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, x, elements[i], E, Poisson, 1.0);
            PlaneMass<double, ShapeFunction8Square, Gauss9Square >(Me, x, elements[i], rho, 1.0);
            Assembling(K, Ke, elements[i], field);
            Assembling(M, Me, elements[i], field);
		}

        //----------Set Dirichlet Boundary Condition----------
        SetDirichlet(K, M, isufixed, ufixed, 1.0e10);

		//----------Solve System Equation----------
        CSR<double> Kmod = CSR<double>(K);
        CSR<double> Mmod = CSR<double>(M);
        std::vector<double> alpha, beta;
        std::vector<std::vector<double> > q;                                                            
        std::vector<double> frequencies = std::vector<double>(m);                                       //Eigenvalues
        std::vector<std::vector<Vector<double> > > u = std::vector<std::vector<Vector<double> > >(m);   //Eigenvectors
        
        LanczosInversePowerProcessForGeneral(Kmod, Mmod, alpha, beta, q, m);
        for(int i = 0; i < m; i++){
            //----------Get eigen value and eigen vector----------
            double eigenvalue = BisectionMethod(alpha, beta, (m - 1) - i);
            frequencies[i] = sqrt(1.0 / eigenvalue);
            std::vector<double> y = InversePowerMethod(alpha, beta, eigenvalue);
            std::vector<double> result = ReconvertVector(y, q);

            //----------Post Process----------
            FieldResultToNodeValue(result, u[i], field);           
        }

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 23), fout);
		AddPointVectors(u[0], "u" + std::to_string(0), fout, true);
		for(int l = 1; l < m; l++){
			AddPointVectors(u[l], "u" + std::to_string(l), fout, false);
		}
		AddElementScalers(s, "s", fout, true);
		fout.close();

		//**************************************************
		//	Get sensitivity and update design variables
		//**************************************************
		double tmp = 0.0;
		for(int l = 0; l < n; l++){
			tmp += 1.0/frequencies[l];
		}
		double objective = (double)n/tmp;											//Function value of objective
		std::vector<double> dobjectives = std::vector<double>(elements.size());		//Sensitivities of objective		
		double weight = 0.0;														//Function values of weight
		std::vector<double> dweights = std::vector<double>(elements.size());		//Sensitivities of weight
		
		//----------Get function values and sensitivities----------
		std::vector<std::vector<double> > dfrequencies = std::vector<std::vector<double> >(elements.size(), std::vector<double>(n));	//Sensitivities of frequency
		std::vector<double> uMu = std::vector<double>(n, 0.0);
		for (int i = 0; i < elements.size(); i++) {
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square>(Ke, x, elements[i], 1.0, Poisson, 1.0);
			Matrix<double> Me;
			PlaneMass<double, ShapeFunction8Square, Gauss9Square>(Me, x, elements[i], 1.0, 1.0);
			
			for(int l = 0; l < n; l++){
				Vector<double> ue = Vector<double>();
				for(int j = 0; j < elements[i].size(); j++){
					ue = ue.Vstack(u[l][elements[i][j]]);
				}
				double ueKeue = ue*(Ke*ue);
				double ueMeue = ue*(Me*ue);
				dfrequencies[i][l] = fabs(p*(E1*pow(s[i], p - 1.0) - E0*pow(s[i], p - 1.0))*ueKeue - frequencies[l]*(rho1 - rho0)*ueMeue);
				if(i == 1389 || i == 1390 || i == 1409 || i == 1410){
					uMu[l] += (rho1*s[i] + rho0*(1.0 - s[i]) + 1.25e-6)*ueMeue;
				} else {
					uMu[l] += (rho1*s[i] + rho0*(1.0 - s[i]))*ueMeue;
				}
			}		

			weight += s[i] - weightlimit;
			if(i == 1389 || i == 1390 || i == 1409 || i == 1410){
				weight += 1.25e-6;
			}
			dweights[i] = 1.0;
		}

		for(int i = 0; i < elements.size(); i++){
			dobjectives[i] = 0.0;
			for(int l = 0; l < n; l++){
				dobjectives[i] += dfrequencies[i][l]/(uMu[l]*pow(frequencies[l], 2.0));
			}
			dobjectives[i] *= (double)n/pow(tmp, 2.0);
		}

		//----------Check convergence----------
		if(fabs((objective - objectivebefore) / (objective + objectivebefore)) < objectiveeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}

		std::cout << "objective:\t" << objective << "\t(\t";
		for(int l = 0; l < n; l++){
			std::cout << frequencies[l] << "\t";
		}
		std::cout << ")\tWeight:\t" << weight << "\t";

		//----------Get updated design variables with OC method----------
		double lambda0 = lambdamin, lambda1 = lambdamax, lambda;
		std::vector<double> snext = std::vector<double>(elements.size());
		while((lambda1 - lambda0) / (lambda1 + lambda0) > lambdaeps){
			lambda = 0.5 * (lambda1 + lambda0);

			for (int i = 0; i < elements.size(); i++) {
				snext[i] = pow(dobjectives[i] / (dweights[i] * lambda), iota) * s[i];
				if(snext[i] < std::max(0.0, (1.0 - movelimit)*s[i])) {
					snext[i] = std::max(0.0, (1.0 - movelimit)*s[i]);
				} else if(snext[i] > std::min(1.0, (1.0 + movelimit)*s[i])) {
					snext[i] = std::min(1.0, (1.0 + movelimit)*s[i]);
				}
			}

			double weightnext = 0.0;
			for (int i = 0; i < elements.size(); i++) {
				if(i == 1389 || i == 1390 || i == 1409 || i == 1410){
					weightnext += 1.25e-6;
				}
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
		objectivebefore = objective;
	}
	
	return 0;
}