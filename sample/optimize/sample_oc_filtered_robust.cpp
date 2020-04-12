#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Equation/Geometric.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/NewtonCotesIntegration.h"
#include "../../src/Optimize/Solver/OC.h"
#include "../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


int main() {
    //----------Define design parameters----------
	double E0 = 0.0001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	double sigmatheta = 5.0*M_PI/180.0;
    double Edtheta2 = pow(sigmatheta, 2.0);
    double Edtheta4 = 3.0*pow(sigmatheta, 4.0);
    double alpha = 1.0;

	double weightlimit = 0.3;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 0.5;


    //----------Generate design region----------
	SquareMesh<double> mesh = SquareMesh<double>(100.0, 100.0, 100, 100);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<int> field = mesh.GenerateFields(2);
    int KDEGREE = *(field.end() - 1);
    std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(fabs(_x(1)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> ufixed = std::vector<double>(isufixed.size(), 0.0);
    std::vector<int> isqfixed0 = mesh.GenerateFixedlist(2, { 1 }, [](Vector<double> _x){
        if(fabs(_x(0) - 50.0) < 1.0e-5 && fabs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed0 = std::vector<double>(isqfixed0.size(), -1.0);
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
	SetNeumann(F, isqfixed0, qfixed0);
    std::vector<int> isqfixed1 = mesh.GenerateFixedlist(2, { 0 }, [](Vector<double> _x){
        if(fabs(_x(0) - 50.0) < 1.0e-5 && fabs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed1 = std::vector<double>(isqfixed1.size(), 1.0);
    std::vector<double> dFdmu = std::vector<double>(KDEGREE, 0.0);
	SetNeumann(dFdmu, isqfixed1, qfixed1);
    std::vector<int> isqfixed2 = mesh.GenerateFixedlist(2, { 1 }, [](Vector<double> _x){
        if(fabs(_x(0) - 50.0) < 1.0e-5 && fabs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed2 = std::vector<double>(isqfixed2.size(), 1.0);
    std::vector<double> d2Fdmu2 = std::vector<double>(KDEGREE, 0.0);
	SetNeumann(d2Fdmu2, isqfixed2, qfixed2);
    

	//----------Get cg of element----------
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        cg[i] = CenterOfGravity(nodes, elements[i]);
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
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);
    OC<double> optimizer = OC<double>(s.size(), 1.0, 1.0e-4, 1.0e4, 1.0e-3, 0.15, std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
			

	//----------Optimize loop----------
	for(int k = 0; k < 500; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%40 == 0){
            beta*=2.0;
            filter.UpdateBeta(beta);
        }


        //----------Get filterd design variables----------
        std::vector<double> rho = filter.GetFilteredVariables(s);


        //----------Get weight value and sensitivities----------
        double g = 0.0;														//Function values of weight
		std::vector<double> dgdrho = std::vector<double>(s.size(), 0.0);    //Sensitivities of weight
        for(int i = 0; i < elements.size(); i++){
            g += scale1*rho[i]/(weightlimit*elements.size());
            dgdrho[i] = scale1/(weightlimit*elements.size()); 
        }
        g -= 1.0*scale1;
        std::vector<double> dgds = filter.GetFilteredSensitivitis(s, dgdrho);

    
        //----------Get compliance value and sensitivities----------
        double f = 0.0;													    //Function value of compliance
		std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);    //Sensitivities of compliance

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		for (int i = 0; i < elements.size(); i++) {
			double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			Matrix<double> Ke;
			PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			Assembling(K, Ke, elements[i], field);
		}

		SetDirichlet(K, isufixed, ufixed, 1.0e10);
        CSR<double> Kmod = CSR<double>(K);	
        std::vector<double> d = ScalingCG(Kmod, F, 100000, 1.0e-10);
        std::vector<double> dddmu = ScalingCG(Kmod, dFdmu, 100000, 1.0e-10);
        std::vector<double> d2ddmu2 = ScalingCG(Kmod, d2Fdmu2, 100000, 1.0e-10);

        double c0 = std::inner_product(F.begin(), F.end(), d.begin(), 0.0);
        double c1 = std::inner_product(dFdmu.begin(), dFdmu.end(), d.begin(), 0.0) + std::inner_product(F.begin(), F.end(), dddmu.begin(), 0.0);
        double c2 = std::inner_product(d2Fdmu2.begin(), d2Fdmu2.end(), d.begin(), 0.0) + 2.0*std::inner_product(dFdmu.begin(), dFdmu.end(), dddmu.begin(), 0.0) + std::inner_product(F.begin(), F.end(), d2ddmu2.begin(), 0.0);

        double EC = c0 + 0.5*c2*Edtheta2;
        double VC = pow(c1, 2.0)*Edtheta2 + 0.25*pow(c2, 2.0)*(Edtheta4 - pow(Edtheta2, 2.0));
        f = scale0*(EC + alpha*sqrt(VC));
        
        double beta1 = alpha*c1*Edtheta2/sqrt(VC);
        double beta2 = 0.5*Edtheta2 + 0.25*alpha*c2*(Edtheta4 - pow(Edtheta2, 2.0))/sqrt(VC);

        std::vector<Vector<double> > dv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d, dv, field);
        std::vector<Vector<double> > dddmuv = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(dddmu, dddmuv, field);
        std::vector<Vector<double> > d2ddmu2v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
		FieldResultToNodeValue(d2ddmu2, d2ddmu2v, field);

        for(int i = 0; i < elements.size(); i++){
            Matrix<double> Ke;
		    PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);

			Vector<double> de = Vector<double>();           
			Vector<double> dddmue = Vector<double>();
			Vector<double> d2ddmu2e = Vector<double>();

            for(int j = 0; j < elements[i].size(); j++){
                de = de.Vstack(dv[elements[i][j]]);
                dddmue = dddmue.Vstack(dddmuv[elements[i][j]]);
                d2ddmu2e = d2ddmu2e.Vstack(d2ddmu2v[elements[i][j]]);
            }

            Vector<double> phi0e = de + beta1*dddmue + beta2*d2ddmu2e;
            Vector<double> phi1e = beta1*de + 2.0*beta2*dddmue;
            Vector<double> phi2e = beta2*de;

            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*((phi0e*(Ke*de)) + (phi1e*(Ke*dddmue)) + (phi2e*(Ke*d2ddmu2e)));
        }

        std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);
        
		       
        //----------Post Process----------
		std::ofstream fout("sample/optimize/result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(dv, "d", fout, true);
		AddElementScalers(rho, "s", fout, true);
		fout.close();
       

		//----------Check convergence----------
        std::cout << "Objective:\t" << f/scale0 << "\tWeight:\t" << g/scale1 << "\t";
		if(optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		

		//----------Get updated design variables with OC----------
		optimizer.UpdateVariables(s, f, dfds, g, dgds, 
            [&](std::vector<double> _xkp1) {
                double g = 0.0;
                std::vector<double> rho = filter.GetFilteredVariables(_xkp1);
                for(int i = 0; i < elements.size(); i++){
                    g += scale1*rho[i]/(weightlimit*elements.size());
                }
                return g - 1.0*scale1; 
            }
        );	
	}
	
	return 0;
}