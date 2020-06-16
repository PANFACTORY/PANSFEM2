#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/FEM/Equation/Truss.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/Optimize/Solver/OC.h"
#include "../../src/PrePost/Mesher/GrandStructure.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    //----------Define design parameters----------
	double E = 210000.0;
	double A0 = 1.0e-5;
    double A1 = 1.5e2;
	double weightlimit = 0.5;
    double scale0 = 1.0e-3;


    //----------Generate design region----------
	GrandStructure2D<double> mesh = GrandStructure2D<double>(4000.0, 2000.0, 4, 2, 1500);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5 && (abs(_x(1)) < 1.0e-5 || abs(_x(1) - 2000.0) < 1.0e-5)) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 4000.0) < 1.0e-5 && abs(_x(1) - 1000.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -5000.0;
    }
    std::vector<double> l = std::vector<double>(elements.size());
    double lsum = 0.0;
    for(int i = 0; i < elements.size(); i++) {
        l[i] = (x[elements[i][1]] - x[elements[i][0]]).Norm();
        lsum += l[i];
    }


	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);
    OC<double> optimizer = OC<double>(s.size(), 0.5, 0.0, 1.0e3, 1.0e-3, 0.15, std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
			

	//----------Optimize loop----------
	for(int k = 0; k < 50; k++){
		std::cout << "\nk = " << k << "\t";
        
        //*************************************************
        //  Get weight value and sensitivities
        //*************************************************
        double g = 0.0;														//Function values of weight
		std::vector<double> dgds = std::vector<double>(s.size(), 0.0);      //Sensitivities of weight
        for(int i = 0; i < elements.size(); i++){
            g += s[i]*l[i]/(weightlimit*lsum);
            dgds[i] = l[i]/(weightlimit*lsum); 
        }
        g -= 1.0;

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        
        //--------------------Get displacement--------------------
		std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        
        SetDirichlet(u, nodetoglobal, ufixed);
        int KDEGREE = Renumbering(nodetoglobal);

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
			double A = A0*(1.0 - s[i]) + A1*s[i];
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Truss2D<double>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, A);
            Assembling(K, F, u, Ke, nodetoglobal, nodetoelement, elements[i]);
		}
        Assembling(F, qfixed, nodetoglobal);

        CSR<double> Kmod = CSR<double>(K);	
        std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(u, result, nodetoglobal);

        //--------------------Get reaction force--------------------
        RemoveBoundaryConditions(nodetoglobal);
        KDEGREE = Renumbering(nodetoglobal);
        std::vector<double> RF = std::vector<double>(KDEGREE, 0.0);
        std::vector<Vector<double> > r = std::vector<Vector<double> >(x.size(), Vector<double>(2));	

        for (int i = 0; i < elements.size(); i++) {
            double A = A0*(1.0 - s[i]) + A1*s[i];
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Truss2D<double>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, A);
            Vector<double> Keue = Ke*ElementVector(u, nodetoelement, elements[i]);
            Assembling(RF, Keue, nodetoglobal, nodetoelement, elements[i]);
        }

        Disassembling(r, RF, nodetoglobal);

        //--------------------Get compliance and sensitivities--------------------
        double f = scale0*std::inner_product(u.begin(), u.end(), r.begin(), 0.0);
		std::vector<double> dfds = std::vector<double>(s.size(), 0.0);
      
        for(int i = 0; i < elements.size(); i++){
            double dA = -A0 + A1;
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Truss2D<double>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, dA);
			Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            dfds[i] = -scale0*ue*(Ke*ue);
        }


        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout("sample/truss/result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 3), fout);
		AddPointVectors(u, "u", fout, true);
        AddPointVectors(r, "r", fout, false);
		AddElementScalers(s, "s", fout, true);
		fout.close();
       

        //*************************************************
        //  Update design variables with OC
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << f/scale0 << "\tWeight:\t" << g << "\t";
		if(optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with OC----------
		optimizer.UpdateVariables(s, f, dfds, g, dgds, 
            [&](std::vector<double> _skp1) {
                double g = 0.0;
                for(int i = 0; i < elements.size(); i++){
                    g += _skp1[i]*l[i]/(weightlimit*lsum);
                }
                return g - 1.0; 
            }
        );	
	}
	
	return 0;
}