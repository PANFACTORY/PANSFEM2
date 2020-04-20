#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/Optimize/Solver/OC.h"
#include "../../src/Optimize/Filter/Heaviside.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    //----------Define design parameters----------
	double E0 = 0.0001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

    double R = 1.5;

    double sigmamu = 5.0*M_PI/180.0;
    double Edmu2 = pow(sigmamu, 2.0);
    double Edmu4 = 3.0*pow(sigmamu, 4.0);
    double alpha = 1.0;

	double weightlimit = 0.3;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 0.5;


    //----------Generate design region----------
	SquareMesh<double> mesh = SquareMesh<double>(100.0, 100.0, 100, 100);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(1)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed0 = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 50.0) < 1.0e-5 && abs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed0) {
        qfixedi.second = -1.0;
    }
    std::vector<std::pair<std::pair<int, int>, double> > qfixed1 = mesh.GenerateFixedlist({ 0 }, [](Vector<double> _x){
        if(abs(_x(0) - 50.0) < 1.0e-5 && abs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed1) {
        qfixedi.second = 1.0;
    }
    std::vector<std::pair<std::pair<int, int>, double> > qfixed2 = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 50.0) < 1.0e-5 && abs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed2) {
        qfixedi.second = 1.0;
    }

	//----------Get cg of element----------
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        cg[i] = CenterOfGravity(x, elements[i]);
    }

    //----------Get neighbor elements list and weight----------
    std::vector<std::vector<int> > neighbors = std::vector<std::vector<int> >(elements.size(), std::vector<int>());
    std::vector<std::vector<double> > w = std::vector<std::vector<double> >(elements.size(), std::vector<double>());
    
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
	std::vector<double> s = std::vector<double>(elements.size(), weightlimit);

	//----------Initialize Optimizer filter class----------
    OC<double> optimizer = OC<double>(s.size(), 1.0, 0.0, 1.0e4, 1.0e-3, 0.15, std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
			
	//----------Optimize loop----------
	for(int k = 0; k < 500; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%40 == 0){
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

        std::vector<double> dgds = filter.GetFilteredSensitivitis(s, dgdrho);

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        
        //--------------------Get displacement--------------------
		std::vector<Vector<double> > u0 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > u1 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > u2 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        
        SetDirichlet(u0, nodetoglobal, ufixed);
        SetDirichlet(u1, nodetoglobal, ufixed);
        SetDirichlet(u2, nodetoglobal, ufixed);
        int KDEGREE = Renumbering(nodetoglobal);

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F0 = std::vector<double>(KDEGREE, 0.0);
        std::vector<double> F1 = std::vector<double>(KDEGREE, 0.0);
        std::vector<double> F2 = std::vector<double>(KDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
			double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, 0.3, 1.0);
            Assembling(K, Ke, nodetoglobal, nodetoelement, elements[i]);
            Assembling(F0, u0, Ke, nodetoglobal, nodetoelement, elements[i]);
            Assembling(F1, u1, Ke, nodetoglobal, nodetoelement, elements[i]);
            Assembling(F2, u2, Ke, nodetoglobal, nodetoelement, elements[i]);
		}
        Assembling(F0, qfixed0, nodetoglobal);
        Assembling(F1, qfixed1, nodetoglobal);
        Assembling(F2, qfixed2, nodetoglobal);

        CSR<double> Kmod = CSR<double>(K);	
        std::vector<double> result0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);
        Disassembling(u0, result0, nodetoglobal);
        std::vector<double> result1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
        Disassembling(u1, result1, nodetoglobal);
        std::vector<double> result2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);
        Disassembling(u2, result2, nodetoglobal);

        //--------------------Get reaction force--------------------
        RemoveBoundaryConditions(nodetoglobal);
        KDEGREE = Renumbering(nodetoglobal);
        F0 = std::vector<double>(KDEGREE, 0.0);
        F1 = std::vector<double>(KDEGREE, 0.0);
        F2 = std::vector<double>(KDEGREE, 0.0);
        std::vector<Vector<double> > r0 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > r1 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > r2 = std::vector<Vector<double> >(x.size(), Vector<double>(2));

        for (int i = 0; i < elements.size(); i++) {
            double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, 0.3, 1.0);
            Vector<double> Keue0 = Ke*ElementVector(u0, nodetoelement, elements[i]);
            Assembling(F0, Keue0, nodetoglobal, nodetoelement, elements[i]);
            Vector<double> Keue1 = Ke*ElementVector(u1, nodetoelement, elements[i]);
            Assembling(F1, Keue1, nodetoglobal, nodetoelement, elements[i]);
            Vector<double> Keue2 = Ke*ElementVector(u2, nodetoelement, elements[i]);
            Assembling(F2, Keue2, nodetoglobal, nodetoelement, elements[i]);
        }

        Disassembling(r0, F0, nodetoglobal);
        Disassembling(r1, F1, nodetoglobal);
        Disassembling(r2, F2, nodetoglobal);

        //--------------------Get compliance and sensitivities--------------------
        double c0 = std::inner_product(r0.begin(), r0.end(), u0.begin(), 0.0);
        double c1 = std::inner_product(r1.begin(), r1.end(), u0.begin(), 0.0) + std::inner_product(r0.begin(), r0.end(), u1.begin(), 0.0);
		double c2 = std::inner_product(r2.begin(), r2.end(), u0.begin(), 0.0) + 2.0*std::inner_product(r1.begin(), r1.end(), u1.begin(), 0.0) + std::inner_product(r0.begin(), r0.end(), u2.begin(), 0.0);
        
        double EC = c0 + 0.5*c2*Edmu2;
        double VC = pow(c1, 2.0)*Edmu2 + 0.25*pow(c2, 2.0)*(Edmu4 - pow(Edmu2, 2.0));
        double f = scale0*(EC + alpha*sqrt(VC));
        
        double beta1 = alpha*c1*Edmu2/sqrt(VC);
        double beta2 = 0.5*Edmu2 + 0.25*alpha*c2*(Edmu4 - pow(Edmu2, 2.0))/sqrt(VC);     
        
        std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);
      
        for(int i = 0; i < elements.size(); i++){
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, 1.0, 0.3, 1.0);
			Vector<double> u0e = ElementVector(u0, nodetoelement, elements[i]);
            Vector<double> u1e = ElementVector(u1, nodetoelement, elements[i]);
            Vector<double> u2e = ElementVector(u2, nodetoelement, elements[i]);

            Vector<double> phi0e = u0e + beta1*u1e + beta2*u2e;
            Vector<double> phi1e = beta1*u0e + 2.0*beta2*u1e;
            Vector<double> phi2e = beta2*u0e;

            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*((phi0e*(Ke*u0e)) + (phi1e*(Ke*u1e)) + (phi2e*(Ke*u2e)));
        }

        std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);


        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout("sample/optimize/result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(u0, "u0", fout, true);
        AddPointVectors(r0, "r0", fout, false);
        AddPointVectors(u1, "u1", fout, false);
        AddPointVectors(r1, "r1", fout, false);
        AddPointVectors(u2, "u2", fout, false);
        AddPointVectors(r2, "r2", fout, false);
		AddElementScalers(rho, "s", fout, true);
		fout.close();
       

        //*************************************************
        //  Update design variables with OC
        //*************************************************

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