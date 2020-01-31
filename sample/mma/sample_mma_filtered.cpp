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
	std::string model_path = "sample/mma/Cantileber3/";
	
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

    //----------Get cg of element----------
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        Vector<double> centeri = Vector<double>(2);
        for(auto nodei : elements[i]){
            centeri += nodes[nodei];
        }
        cg[i] = centeri/(double)elements[i].size();
    }

    //----------Get neighbor elements list and weight----------
    std::vector<std::vector<int> > neighbors = std::vector<std::vector<int> >(elements.size(), std::vector<int>());
    std::vector<std::vector<double> > w = std::vector<std::vector<double> >(elements.size(), std::vector<double>());
    double R = 1.5;
    for(int i = 0; i < elements.size(); i++){
        for(int j = 0; j < elements.size(); j++){
            if((cg[i] - cg[j]).Norm() <= R){
                neighbors[i].push_back(j);
                w[i].push_back((R - (cg[i] - cg[j]).Norm())/R);
            }
        }
    }

	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.0001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	double weightlimit = 0.5;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 1.0;

    MMA<double> optimizer = MMA<double>(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.2, 0.5, 0.7, 1.2, 1.0e-6);
		
	//----------Optimize loop----------
	for(int k = 0; k < 500; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%40 == 1){
            beta*=2.0;
        }

        //*************************************************
        //  Get filterd design variables
        //*************************************************
        std::vector<double> stilde = std::vector<double>(s.size());
        std::vector<double> rho = std::vector<double>(s.size());
        for(int i = 0; i < elements.size(); i++){
            double wssum = 0.0;
            double wsum = 0.0;
            for(int j = 0; j < neighbors[i].size(); j++){
                wssum += w[i][j]*s[neighbors[i][j]];
                wsum += w[i][j];
            }
            stilde[i] = wssum/wsum;
            rho[i] = 0.5*(tanh(0.5*beta) + tanh(beta*(stilde[i] - 0.5)))/tanh(0.5*beta);
        }


        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        double objective = 0.0;													//Function value of compliance
		std::vector<double> dobjectives = std::vector<double>(s.size(), 0.0);	//Sensitivities of compliance

        std::vector<double> constraints = std::vector<double>(1);																//Function values of weight
		std::vector<std::vector<double> > dconstraints = std::vector<std::vector<double> >(1, std::vector<double>(s.size(), 0.0));	//Sensitivities of weight

        //----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1 * pow(rho[i], p) + E0 * (1.0 - pow(rho[i], p));
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
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
        
        std::vector<double> dfdrho = std::vector<double>(s.size());
        std::vector<double> drhodstilde = std::vector<double>(s.size());
        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);
			Vector<double> de = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                de = de.Vstack(dv[elements[i][j]]);
            }
            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*(de*(Ke*de));
            drhodstilde[i] = 0.5*beta*(1.0 - pow(tanh(beta*(stilde[i] - 0.5)), 2.0))/tanh(0.5*beta);

            constraints[0] += scale1*rho[i]/(weightlimit*elements.size());
        }
        constraints[0] -= 1.0*scale1;

        for(int i = 0; i < elements.size(); i++){
            double wsum = 0.0;
            for(int j = 0; j < neighbors[i].size(); j++){
                dobjectives[i] += dfdrho[neighbors[i][j]]*drhodstilde[neighbors[i][j]]*w[i][j];
                dconstraints[0][i] += scale1*drhodstilde[neighbors[i][j]]*w[i][j]/(weightlimit*elements.size());
                wsum += w[i][j];
            }
            dobjectives[i] /= wsum;      
            dconstraints[0][i] /= wsum; 
        }
		
        
        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(dv, "d", fout, true);
		AddElementScalers(rho, "s", fout, true);
		fout.close();
       

        //*************************************************
        //  Update design variables with OC method
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << objective/scale0 << "\tWeight:\t" << constraints[0]/scale1 << "\t";
		if(optimizer.IsConvergence(objective)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with OC method----------
		optimizer.UpdateVariables(s, objective, dobjectives, constraints, dconstraints);	
	}
	
	return 0;
}