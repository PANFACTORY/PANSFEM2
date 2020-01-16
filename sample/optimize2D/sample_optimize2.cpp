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
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/optimize2D/";
	
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
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
	SetNeumann(F, isqfixed, qfixed);
	
	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(2*elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 823.0;
    double E2 = 210000.0;
	double Poisson = 0.3;
    double rho0 = 0.0;
    double rho1 = 0.0323;
    double rho2 = 1.0;

	double p = 3.0;
    double q = 3.0;

	double iota = 0.75;
	double lambdamin = 1.0e-20;
	double lambdamax = 1.0e20;
	double lambdaeps = 1.0e-15;
	double movelimit = 0.15;

	double weightlimit = 0.05;
	double objectivebefore = 0.0;
	double objectiveeps = 1.0e-5;
	
	//----------Optimize loop----------
	for(int k = 0; k < 200; k++){
		std::cout << "k = " << k << "\t";

        
        //*************************************************
		//  Get weight value and sensitivities
		//*************************************************
        double weight = 0.0;													//Function values of weight
		std::vector<double> dweights = std::vector<double>(s.size());			//Sensitivities of weight
        for (int i = 0; i < elements.size(); i++) {						
			weight += rho0*(1.0 - s[2*i]) + (rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])*s[2*i] - weightlimit;
			dweights[2*i] = - rho0 + rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1];
            dweights[2*i + 1] = (- rho1 + rho2)*s[2*i];
		}

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        double objective = 0.0;													//Function value of compliance
		std::vector<double> dobjectives = std::vector<double>(s.size(), 0.0);	//Sensitivities of compliance

        //----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E0*(1.0 - pow(s[2*i], p)) + (E1*(1.0 - pow(s[2*i + 1], q)) + E2*pow(s[2*i + 1], q))*pow(s[2*i], p);
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			Assembling(K, Ke, elements[i], field);
		}

		//----------Set Dirichlet Boundary Condition(Only fix displacement 0.0)----------
		SetDirichlet(K, isufixed, ufixed, 1.0e10);
        CSR<double> Kmod = CSR<double>(K);	

        //----------Get function value and sensitivities----------
        std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);
        
        objective = std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
       
        std::vector<Vector<double> > dv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d, dv, field);

        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);

			Vector<double> de = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                de = de.Vstack(dv[elements[i][j]]);
            }

            dobjectives[2*i] = p*(-E0 + (E1*(1.0 - pow(s[2*i + 1], q)) + E2*pow(s[2*i + 1], q)))*pow(s[2*i], p - 1.0)*(de*(Ke*de));
            dobjectives[2*i + 1] = q*(-E1 + E2)*pow(s[2*i + 1], q - 1.0)*pow(s[2*i], p)*(de*(Ke*de));
        }

        
        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 23), fout);
		AddPointVectors(dv, "d", fout, true);
		std::vector<double> s0 = std::vector<double>(elements.size());
        std::vector<double> s1 = std::vector<double>(elements.size());
        std::vector<double> rho = std::vector<double>(elements.size());
        for(int i = 0; i < elements.size(); i++){
            s0[i] = s[2*i];
            s1[i] = s[2*i + 1];
            rho[i] = rho0*(1.0 - s[2*i]) + (rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])*s[2*i];
        }
		AddElementScalers(s0, "s0", fout, true);
        AddElementScalers(s1, "s1", fout, false);
        AddElementScalers(rho, "rho", fout, false);
		fout.close();
       

        //*************************************************
        //  Update design variables with OC method
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << objective << "\tWeight:\t" << weight << "\t";
		if(fabs((objective - objectivebefore) / (objective + objectivebefore)) < objectiveeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with OC method----------
		double lambda0 = lambdamin, lambda1 = lambdamax, lambda;
		std::vector<double> snext = std::vector<double>(s.size());
		while((lambda1 - lambda0) / (lambda1 + lambda0) > lambdaeps){
			lambda = 0.5 * (lambda1 + lambda0);

			for (int i = 0; i < s.size(); i++) {
				snext[i] = pow(dobjectives[i] / (dweights[i] * lambda), iota) * s[i];
				if(snext[i] < std::max(0.0, (1.0 - movelimit)*s[i])) {
					snext[i] = std::max(0.0, (1.0 - movelimit)*s[i]);
				} else if(snext[i] > std::min(1.0, (1.0 + movelimit)*s[i])) {
					snext[i] = std::min(1.0, (1.0 + movelimit)*s[i]);
				}
			}

			double weightnext = 0.0;
			for (int i = 0; i < elements.size(); i++) {
				weightnext += rho0*(1.0 - snext[2*i]) + (rho1*(1.0 - snext[2*i + 1]) + rho2*snext[2*i + 1])*snext[2*i] - weightlimit;
			}

			if (weightnext > 0.0) {
				lambda0 = lambda;
			}
			else {
				lambda1 = lambda;
			}
		}

		std::cout << "Lagrange value:\t" << lambda << std::endl;

		//----------Update design variables and objective----------
		for (int i = 0; i < s.size(); i++) {
			s[i] = snext[i];
		}
		objectivebefore = objective;
	}
	
	return 0;
}