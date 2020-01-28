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

    std::vector<int> isqfixed2;
	std::vector<double> qfixed2;
	ImportNeumannFromCSV(isqfixed2, qfixed2, field, model_path + "Neumann2.csv");
    std::vector<double> dFdtheta = std::vector<double>(KDEGREE, 0.0);
    SetNeumann(dFdtheta, isqfixed2, qfixed2);

    std::vector<int> isqfixed3;
	std::vector<double> qfixed3;
	ImportNeumannFromCSV(isqfixed3, qfixed3, field, model_path + "Neumann3.csv");
    std::vector<double> d2Fdtheta2 = std::vector<double>(KDEGREE, 0.0);
    SetNeumann(d2Fdtheta2, isqfixed3, qfixed3);
	
	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(2*elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 23.3;
    double E2 = 210000.0;
	double Poisson = 0.3;
    double rho0 = 0.0;
    double rho1 = 0.0323;
    double rho2 = 1.0;

	double p = 3.0;
    double q = 3.0;

    double sigmatheta = 15.0/180.0*3.141592;
    double Edtheta2 = pow(sigmatheta, 2.0);
    double Edtheta4 = 3.0*pow(sigmatheta, 4.0);
    double alpha = 1.0;

	double weightlimit = 0.2;
    double voidlimit = 0.3;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

	MMA<double> optimizer = MMA<double>(s.size(), 2, 1.0,
		std::vector<double>(2, 0.0),
		std::vector<double>(2, 10000.0),
		std::vector<double>(2, 0.0), 
		std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.2, 0.5, 0.7, 1.2, 1.0e-5);
			
	//----------Optimize loop----------
	for(int k = 0; k < 200; k++){
		std::cout << "\nk = " << k << "\t";

        
        //*************************************************
		//  Get weight value and sensitivities
		//*************************************************
        std::vector<double> constraints = std::vector<double>(2);	//Function values of weight
		std::vector<std::vector<double> > dconstraints = std::vector<std::vector<double> >(2, std::vector<double>(s.size()));	//Sensitivities of weight
        for (int i = 0; i < elements.size(); i++) {						
			constraints[0] += (rho0*(1.0 - s[2*i]) + (rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])*s[2*i])/(weightlimit*elements.size());
			dconstraints[0][2*i] = (- rho0 + rho1*(1.0 - s[2*i + 1]) + rho2*s[2*i + 1])/(weightlimit*elements.size());
            dconstraints[0][2*i + 1] = (- rho1 + rho2)*s[2*i]/(weightlimit*elements.size());

            constraints[1] += (1.0 - s[2*i])/(voidlimit*elements.size());
            dconstraints[1][2*i] = -1.0/(voidlimit*elements.size());
            dconstraints[1][2*i + 1] = 0.0;
		}
		constraints[0] -= 1.0;
        constraints[1] -= 1.0;
        
        //*************************************************
        //  Get robust compliance value and sensitivities
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
        std::vector<double> dddtheta = ScalingCG(Kmod, dFdtheta, 100000, 1.0e-10);
        std::vector<double> d2ddtheta2 = ScalingCG(Kmod, d2Fdtheta2, 100000, 1.0e-10);

        double c0 = std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
        double c1 = std::inner_product(dFdtheta.begin(), dFdtheta.end(), d.begin(), 0.0) + std::inner_product(F.begin(), F.end(), dddtheta.begin(), 0.0);
        double c2 = std::inner_product(d2Fdtheta2.begin(), d2Fdtheta2.end(), d.begin(), 0.0) + 2.0*std::inner_product(dFdtheta.begin(), dFdtheta.end(), dddtheta.begin(), 0.0) + std::inner_product(F.begin(), F.end(), d2ddtheta2.begin(), 0.0);

        double EC = c0 + 0.5*c2*Edtheta2;
        double VC = pow(c1, 2.0)*Edtheta2 + 0.25*pow(c2, 2.0)*(Edtheta4 - pow(Edtheta2, 2.0));
        objective = scale0*(EC + alpha*sqrt(VC));
        
        double beta1 = alpha*c1*Edtheta2/sqrt(VC);
        double beta2 = 0.5*Edtheta2 + 0.25*alpha*c2*(Edtheta4 - pow(Edtheta2, 2.0))/sqrt(VC);

        std::vector<double> F0 = std::vector<double>(KDEGREE);
        std::vector<double> F1 = std::vector<double>(KDEGREE);
        std::vector<double> F2 = std::vector<double>(KDEGREE);
        for(int i = 0; i < KDEGREE; i++){
            F2[i] = beta2*F[i];
            F1[i] = beta1*F[i] + 2.0*beta2*dFdtheta[i];
            F0[i] = F[i] + beta1*dFdtheta[i] + beta2*d2Fdtheta2[i];
        }
        std::vector<double> phi2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);
        std::vector<double> phi1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
        std::vector<double> phi0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);

        std::vector<Vector<double> > phi0v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(phi0, phi0v, field);
        std::vector<Vector<double> > dv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d, dv, field);
        std::vector<Vector<double> > phi1v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(phi1, phi1v, field);
        std::vector<Vector<double> > dddthetav = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(dddtheta, dddthetav, field);
        std::vector<Vector<double> > phi2v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(phi2, phi2v, field);
        std::vector<Vector<double> > d2ddtheta2v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d2ddtheta2, d2ddtheta2v, field);

        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);

            Vector<double> phi0e = Vector<double>();
			Vector<double> de = Vector<double>();
            Vector<double> phi1e = Vector<double>();
			Vector<double> dddthetae = Vector<double>();
            Vector<double> phi2e = Vector<double>();
			Vector<double> d2ddtheta2e = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                phi0e = phi0e.Vstack(phi0v[elements[i][j]]);
                de = de.Vstack(dv[elements[i][j]]);
                phi1e = phi1e.Vstack(phi1v[elements[i][j]]);
                dddthetae = dddthetae.Vstack(dddthetav[elements[i][j]]);
                phi2e = phi2e.Vstack(phi2v[elements[i][j]]);
                d2ddtheta2e = d2ddtheta2e.Vstack(d2ddtheta2v[elements[i][j]]);
            }

            dobjectives[2*i] = p*scale0*(-E0 + (E1*(1.0 - pow(s[2*i + 1], q)) + E2*pow(s[2*i + 1], q)))*pow(s[2*i], p - 1.0)*(- (phi0e*(Ke*de)) - (phi1e*(Ke*dddthetae)) - (phi2e*(Ke*d2ddtheta2e)));
            dobjectives[2*i + 1] = q*scale0*(-E1 + E2)*pow(s[2*i + 1], q - 1.0)*pow(s[2*i], p)*(- (phi0e*(Ke*de)) - (phi1e*(Ke*dddthetae)) - (phi2e*(Ke*d2ddtheta2e)));
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