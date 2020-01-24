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
#include "../../src/Optimize/Solver/MMA.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/mma/";
	
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
	std::vector<double> s = std::vector<double>(elements.size()*2);
	for(int i = 0; i < elements.size(); i++){
		s[2*i] = 0.5;
		s[2*i + 1] = 0.5;
	}
	
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

	double weightlimit = 0.5;
	double scale0 = 1.0e2;
	double scale1 = 1.0;
	
    MMA<double> optimizer = MMA<double>(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.05, 0.5, 0.7, 1.2, 1.0e-5);
		
	//----------Optimize loop----------
	for(int k = 0; k < 100; k++){
		std::cout << "\nk = " << k << "\t";

        
        //*************************************************
		//  Get weight value and sensitivities
		//*************************************************
        std::vector<double> constraints = std::vector<double>(1);																//Function values of weight
		std::vector<std::vector<double> > dconstraints = std::vector<std::vector<double> >(1, std::vector<double>(s.size()));	//Sensitivities of weight
        for (int i = 0; i < elements.size(); i++) {						
			constraints[0] += scale1*(rho0*(1.0 - s[2*i]) + (rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])*s[2*i])/(weightlimit*elements.size());
            dconstraints[0][2*i] = scale1*(- rho0 + rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])/(weightlimit*elements.size());
			dconstraints[0][2*i + 1] = scale1*(- rho1 + rho2)*s[2*i]/(weightlimit*elements.size());
		}
		constraints[0] -= 1.0*scale1;

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        double objective = 0.0;													//Function value of compliance
		std::vector<double> dobjectives = std::vector<double>(s.size(), 0.0);	//Sensitivities of compliance

        //----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1 * pow(s[i], p) + E0 * (1.0 - pow(s[i], p));
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			Assembling(K, Ke, elements[i], field);
		}

		//----------Set Dirichlet Boundary Condition(Only fix displacement 0.0)----------
		SetDirichlet(K, isufixed, ufixed, 1.0e10);
        CSR<double> Kmod = CSR<double>(K);	

        //----------Get function value and sensitivities----------
        std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);
        
        objective = scale0*std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
        
        std::vector<Vector<double> > dv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d, dv, field);
        
        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);

			Vector<double> de = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                de = de.Vstack(dv[elements[i][j]]);
            }

			dobjectives[2*i] = -scale0*p*(-E0 + (E1*(1.0 - pow(s[2*i + 1], q)) + E2*pow(s[2*i + 1], q)))*pow(s[2*i], p - 1.0)*(de*(Ke*de));
            dobjectives[2*i + 1] = -scale0*q*(-E1 + E2)*pow(s[2*i + 1], q - 1.0)*pow(s[2*i], p)*(de*(Ke*de));
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

		std::vector<double> sensc0 = std::vector<double>(elements.size());
		std::vector<double> sensc1 = std::vector<double>(elements.size());
		std::vector<double> sensv0 = std::vector<double>(elements.size());
		std::vector<double> sensv1 = std::vector<double>(elements.size());

        for(int i = 0; i < elements.size(); i++){
            s0[i] = s[2*i];
            s1[i] = s[2*i + 1];
            rho[i] = rho0*(1.0 - s[2*i]) + (rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])*s[2*i];

			sensc0[i] = dobjectives[2*i];
			sensc1[i] = dobjectives[2*i + 1];
			sensv0[i] = dconstraints[0][2*i];
			sensv1[i] = dconstraints[0][2*i + 1];
        }
		AddElementScalers(s0, "s0", fout, true);
        AddElementScalers(s1, "s1", fout, false);
        AddElementScalers(rho, "rho", fout, false);

		AddElementScalers(sensc0, "sensc0", fout, false);
		AddElementScalers(sensc1, "sensc1", fout, false);
		AddElementScalers(sensv0, "sensv0", fout, false);
		AddElementScalers(sensv1, "sensv1", fout, false);
		fout.close();
       

        //*************************************************
        //  Update design variables with OC method
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << objective/scale0 << "\tWeight:\t" << constraints[0]/scale1 << "\t";
		if(optimizer.IsConvergence(objective)){
			//std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			//break;
		}
		
		//----------Get updated design variables with OC method----------
		optimizer.UpdateVariables(s, objective, dobjectives, constraints, dconstraints);	
	}
	
	return 0;
}