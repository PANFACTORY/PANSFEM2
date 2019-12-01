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

	double iota = 0.75;
	double lambdamin = 1.0e-20;
	double lambdamax = 1.0e20;
	double lambdaeps = 1.0e-15;
	double movelimit = 0.1;

	double weightlimit = 0.36;
	double frequencybefore = 0.0;
	double frequencyeps = 1.0e-5;
	
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
			if(i == 3974 || i == 3975){
				rho += 2.0e-6;
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
        for(int i = 0; i < 1; i++){
            //----------Get eigen value and eigen vector----------
            double eigenvalue = BisectionMethod(alpha, beta, (m - 1) - i);
            frequencies[i] = sqrt(1.0 / eigenvalue);
            std::vector<double> y = InversePowerMethod(alpha, beta, eigenvalue);
            std::vector<double> result = ReconvertVector(y, q);

            //----------Post Process----------
            FieldResultToNodeValue(result, u[i], field);

            //----------Save file----------
            std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
            MakeHeadderToVTK(fout);
            AddPointsToVTK(x, fout);
            AddElementToVTK(elements, fout);
            AddElementTypes(std::vector<int>(elements.size(), 23), fout);
            AddPointVectors(u[i], "u", fout);
			AddElementScalers(s, "s", fout);
            fout.close();
        }

		//**************************************************
		//	Get sensitivity and update design variables
		//**************************************************
		double frequency = frequencies[0];											//Function value of frequency
		std::vector<double> dfrequencies = std::vector<double>(elements.size());	//Sensitivities of frequency
		double weight = 0.0;														//Function values of weight
		std::vector<double> dweights = std::vector<double>(elements.size());		//Sensitivities of weight
		double uMu = 0.0;

		//----------Get function values and sensitivities----------
		for (int i = 0; i < elements.size(); i++) {
			Vector<double> ue = Vector<double>();
			for(int j = 0; j < elements[i].size(); j++){
				ue = ue.Vstack(u[0][elements[i][j]]);
			}
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square>(Ke, x, elements[i], 1.0, Poisson, 1.0);
			double ueKeue = ue*(Ke*ue);

			Matrix<double> Me;
			PlaneMass<double, ShapeFunction8Square, Gauss9Square>(Me, x, elements[i], 1.0, 1.0);
			double ueMeue = ue*(Me*ue);
			dfrequencies[i] = fabs(p*(E1*pow(s[i], p - 1.0) - E0*pow(s[i], p - 1.0))*ueKeue - frequencies[0]*(rho1 - rho0)*ueMeue);
			weight += s[i] - weightlimit;
			dweights[i] = 1.0;
			if(i == 3974 || i == 3975){
				uMu += (rho1*s[i] + rho0*(1.0 - s[i]) + 2.0e-6)*ueMeue;
			} else {
				uMu += (rho1*s[i] + rho0*(1.0 - s[i]))*ueMeue;
			}
		}

		for(int i = 0; i < elements.size(); i++){
			dfrequencies[i] /= uMu;
		}

		//----------Check convergence----------
		if(fabs((frequency - frequencybefore) / (frequency + frequencybefore)) < frequencyeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}

		std::cout << "frequency:\t" << frequency << "\t";
		std::cout << "Weight:\t" << weight << "\t";

		//----------Get updated design variables with OC method----------
		double lambda0 = lambdamin, lambda1 = lambdamax, lambda;
		std::vector<double> snext = std::vector<double>(elements.size());
		while((lambda1 - lambda0) / (lambda1 + lambda0) > lambdaeps){
			lambda = 0.5 * (lambda1 + lambda0);

			for (int i = 0; i < elements.size(); i++) {
				snext[i] = pow(dfrequencies[i] / (dweights[i] * lambda), iota) * s[i];
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
		frequencybefore = frequency;
	}
	
	return 0;
}