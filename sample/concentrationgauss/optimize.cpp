#define _USE_MATH_DEFINES
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
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/NewtonCotesIntegration.h"
#include "../../src/Optimize/Solver/MMA.h"
#include "../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


Vector<double> P(Vector<double> _x){
	double sigma = 10.0;
	double mu = 300.0;
    return { 0.0, -exp(-0.5*pow((_x(0) - mu)/sigma, 2.0))/(sqrt(2.0*M_PI)*sigma) };
}


Vector<double> dPdmu(Vector<double> _x){
	double sigma = 10.0;
	double mu = 300.0;
    return { 0.0, -(_x(0) - mu)*exp(-0.5*pow((_x(0) - mu)/sigma, 2.0))/(sqrt(2.0*M_PI)*pow(sigma, 3.0)) };
}


Vector<double> d2Pdmu2(Vector<double> _x){
	double sigma = 10.0;
	double mu = 300.0;
    return { 0.0, -(pow((_x(0) - mu)/sigma, 2.0) - 1.0)*exp(-0.5*pow((_x(0) - mu)/sigma, 2.0))/(sqrt(2.0*M_PI)*pow(sigma, 3.0)) };
}


int main() {
	//----------Model Path----------
	std::string model_path = "sample/concentrationgauss/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

    //----------Add Edge----------
	std::vector<std::vector<int> > edges;
	ImportElementsFromCSV(edges, model_path + "Edge.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");
	
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
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);

	//----------Define design parameters----------
	double E0 = 0.001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	double sigmamu = 100.0;
    double Edmu2 = pow(sigmamu, 2.0);
    double Edmu4 = 3.0*pow(sigmamu, 4.0);
    double alpha = 1.0;

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
	for(int k = 0; k < 200; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%40 == 1){
            beta*=2.0;
            filter.UpdateBeta(beta);
        }

        
        //*************************************************
        //  Get filterd design variables
        //*************************************************
        std::vector<double> rho = filter.GetFilteredVariables(s);


        //*************************************************
        //  Get weight value and sensitivities
        //*************************************************
        double g = 0.0;														//Function values of weight
		std::vector<double> dgdrho = std::vector<double>(s.size(), 0.0);    //Sensitivities of weight
        for(int i = 0; i < elements.size(); i++){
            g += scale1*rho[i]/(weightlimit*elements.size());
            dgdrho[i] = scale1/(weightlimit*elements.size()); 
        }
        g -= 1.0*scale1;

        
        //*************************************************
        //  Get robust compliance value and sensitivities
        //*************************************************
        double f = 0.0;													//Function value of compliance
		std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);	//Sensitivities of compliance

        //----------Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			Matrix<double> Ke;
			PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			Assembling(K, Ke, elements[i], field);
		}

		//----------Set Dirichlet Boundary Condition(Only fix displacement 0.0)----------
		SetDirichlet(K, isufixed, ufixed, 1.0e10);
        CSR<double> Kmod = CSR<double>(K);	

        //----------Get function value and sensitivities----------
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
        for(auto edge : edges){
            Vector<double> Fe;
            PlaneSurfaceForce<double, ShapeFunction2Line, NewtonCotes3Line>(Fe, nodes, edge, P, 1.0);
            Assembling(F, Fe, edge, field);
        }
        std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);

        std::vector<double> dFdmu = std::vector<double>(KDEGREE, 0.0);		//External load vector
        for(auto edge : edges){
            Vector<double> Fe;
            PlaneSurfaceForce<double, ShapeFunction2Line, NewtonCotes3Line>(Fe, nodes, edge, dPdmu, 1.0);
            Assembling(dFdmu, Fe, edge, field);
        }
        std::vector<double> dddmu = ScalingCG(Kmod, dFdmu, 100000, 1.0e-10);

        std::vector<double> d2Fdmu2 = std::vector<double>(KDEGREE, 0.0);		//External load vector
        for(auto edge : edges){
            Vector<double> Fe;
            PlaneSurfaceForce<double, ShapeFunction2Line, NewtonCotes3Line>(Fe, nodes, edge, d2Pdmu2, 1.0);
            Assembling(d2Fdmu2, Fe, edge, field);
        }
        std::vector<double> d2ddmu2 = ScalingCG(Kmod, d2Fdmu2, 100000, 1.0e-10);

        double c0 = std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
        double c1 = std::inner_product(dFdmu.begin(), dFdmu.end(), d.begin(), 0.0) + std::inner_product(F.begin(), F.end(), dddmu.begin(), 0.0);
        double c2 = std::inner_product(d2Fdmu2.begin(), d2Fdmu2.end(), d.begin(), 0.0) + 2.0*std::inner_product(dFdmu.begin(), dFdmu.end(), dddmu.begin(), 0.0) + std::inner_product(F.begin(), F.end(), d2ddmu2.begin(), 0.0);

        double EC = c0 + 0.5*c2*Edmu2;
        double VC = pow(c1, 2.0)*Edmu2 + 0.25*pow(c2, 2.0)*(Edmu4 - pow(Edmu2, 2.0));
        f = scale0*(EC + alpha*sqrt(VC));
        
        double beta1 = alpha*c1*Edmu2/sqrt(VC);
        double beta2 = 0.5*Edmu2 + 0.25*alpha*c2*(Edmu4 - pow(Edmu2, 2.0))/sqrt(VC);

        std::vector<double> F0 = std::vector<double>(KDEGREE);
        std::vector<double> F1 = std::vector<double>(KDEGREE);
        std::vector<double> F2 = std::vector<double>(KDEGREE);
        for(int i = 0; i < KDEGREE; i++){
            F2[i] = beta2*F[i];
            F1[i] = beta1*F[i] + 2.0*beta2*dFdmu[i];
            F0[i] = F[i] + beta1*dFdmu[i] + beta2*d2Fdmu2[i];
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
        std::vector<Vector<double> > dddmuv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(dddmu, dddmuv, field);
        std::vector<Vector<double> > phi2v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(phi2, phi2v, field);
        std::vector<Vector<double> > d2ddmu2v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d2ddmu2, d2ddmu2v, field);

        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrain<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);

            Vector<double> phi0e = Vector<double>();
			Vector<double> de = Vector<double>();
            Vector<double> phi1e = Vector<double>();
			Vector<double> dddmue = Vector<double>();
            Vector<double> phi2e = Vector<double>();
			Vector<double> d2ddmu2e = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                phi0e = phi0e.Vstack(phi0v[elements[i][j]]);
                de = de.Vstack(dv[elements[i][j]]);
                phi1e = phi1e.Vstack(phi1v[elements[i][j]]);
                dddmue = dddmue.Vstack(dddmuv[elements[i][j]]);
                phi2e = phi2e.Vstack(phi2v[elements[i][j]]);
                d2ddmu2e = d2ddmu2e.Vstack(d2ddmu2v[elements[i][j]]);
            }

            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*((phi0e*(Ke*de)) + (phi1e*(Ke*dddmue)) + (phi2e*(Ke*d2ddmu2e)));
        }


        //*************************************************
        //  Filtering sensitivities
        //*************************************************
        std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);
        std::vector<double> dgds = filter.GetFilteredSensitivitis(s, dgdrho);
		
        
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
        //  Update design variables with MMA
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << f/scale0 << "\t(" << EC << "\t" << VC << ")\tWeight:\t" << g/scale1 << "\t";
		if(optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with MMA----------
		optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });	
	}
	
	return 0;
}