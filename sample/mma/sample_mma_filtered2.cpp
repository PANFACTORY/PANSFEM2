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
#include "../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/mma/QBBBeam2/";
	
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

    //----------Initialize Heaviside filter class----------
    HeavisideFilter<double> filter = HeavisideFilter<double>(elements.size(), neighbors, w);

	//----------Initialize design variables----------
	std::vector<double> s0 = std::vector<double>(elements.size(), 0.5);
    std::vector<double> s1 = std::vector<double>(elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 3417.0;
    double E2 = 210000.0;
	double Poisson = 0.3;
    double rho0 = 0.0;
    double rho1 = 0.119;
    double rho2 = 1.0;
	double p = 3.0;
    double q = 3.0;

	double weightlimit = 0.3;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 1.0;

    MMA<double> optimizer = MMA<double>(s0.size() + s1.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s0.size() + s1.size(), 0.01), std::vector<double>(s0.size() + s1.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.2, 0.5, 0.7, 1.2, 1.0e-6);
		
	//----------Optimize loop----------
	for(int k = 0; k < 500; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%40 == 1){
            beta*=2.0;
            filter.UpdateBeta(beta);
        }

        //*************************************************
        //  Get filterd design variables
        //*************************************************
        std::vector<double> r0 = filter.GetFilteredVariables(s0);
        std::vector<double> r1 = filter.GetFilteredVariables(s1);


        //*************************************************
        //  Get weight value and sensitivities
        //*************************************************
        double g = 0.0;														//Function values of weight
		std::vector<double> dgdr0 = std::vector<double>(s0.size(), 0.0);    //Sensitivities of weight
        std::vector<double> dgdr1 = std::vector<double>(s1.size(), 0.0);    //Sensitivities of weight
        for (int i = 0; i < elements.size(); i++) {						
			g += (rho0*(1.0 - r0[i]) + (rho1*(1.0 - r1[i]) + rho2*r1[i])*r0[i])/(weightlimit*elements.size());
			dgdr0[i] = (- rho0 + rho1*(1.0 - r1[i]) + rho2*r1[i])/(weightlimit*elements.size());
            dgdr1[i] = (- rho1 + rho2)*r0[i]/(weightlimit*elements.size());
		}
		g -= 1.0;

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        double f = 0.0;													    //Function value of compliance
		std::vector<double> dfdr0 = std::vector<double>(s0.size(), 0.0);    //Sensitivities of compliance
        std::vector<double> dfdr1 = std::vector<double>(s1.size(), 0.0);    //Sensitivities of compliance

        //----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E0*(1.0 - pow(r0[i], p)) + (E1*(1.0 - pow(r1[i], q)) + E2*pow(r1[i], q))*pow(r0[i], p);
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			Assembling(K, Ke, elements[i], field);
		}

		//----------Set Dirichlet Boundary Condition(Only fix displacement 0.0)----------
		SetDirichlet(K, isufixed, ufixed, 1.0e10);
        CSR<double> Kmod = CSR<double>(K);	

        //----------Get function value and sensitivities----------
        std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);

        f = scale0*std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
        
        std::vector<Vector<double> > dv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d, dv, field);
        
        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);
			Vector<double> de = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                de = de.Vstack(dv[elements[i][j]]);
            }
            dfdr0[i] = -scale0*p*(-E0 + (E1*(1.0 - pow(r1[i], q)) + E2*pow(r1[i], q)))*pow(r0[i], p - 1.0)*(de*(Ke*de));
            dfdr1[i] = -scale0*q*(-E1 + E2)*pow(r1[i], q - 1.0)*pow(r0[i], p)*(de*(Ke*de));
        }


        //*************************************************
        //  Filtering sensitivities
        //*************************************************
        std::vector<double> dfds0 = filter.GetFilteredSensitivitis(s0, dfdr0);
        std::vector<double> dfds1 = filter.GetFilteredSensitivitis(s1, dfdr1);
        std::vector<double> dgds0 = filter.GetFilteredSensitivitis(s0, dgdr0);
        std::vector<double> dgds1 = filter.GetFilteredSensitivitis(s1, dgdr1);
		
        
        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(dv, "d", fout, true);
		AddElementScalers(r0, "r0", fout, true);
        AddElementScalers(r1, "r1", fout, false);
        std::vector<double> rho = std::vector<double>(elements.size());
        for(int i = 0; i < elements.size(); i++){
            rho[i] = rho0*(1.0 - r0[i]) + (rho1*(1.0 - r1[i]) + rho2*r1[i])*r0[i];
        }
        AddElementScalers(rho, "rho", fout, false);
		fout.close();
       

        //*************************************************
        //  Update design variables with MMA
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << f/scale0 << "\tWeight:\t" << g/scale1 << "\t";
		if(optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with MMA----------
        std::vector<double> s = std::vector<double>(s0.size() + s1.size());
        std::vector<double> dfds = std::vector<double>(s0.size() + s1.size());
        std::vector<double> dgds = std::vector<double>(s0.size() + s1.size());
        for(int i = 0; i < elements.size(); i++){
            s[2*i] = s0[i];
            s[2*i + 1] = s1[i];
            dfds[2*i] = dfds0[i];
            dfds[2*i + 1] = dfds1[i];
            dgds[2*i] = dgds0[i];
            dgds[2*i + 1] = dgds1[i];
        }
		optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });	
        for(int i = 0; i < elements.size(); i++){
            s0[i] = s[2*i];
            s1[i] = s[2*i + 1];
        }
	}
	
	return 0;
}