#include <iostream>
#include <vector>


#include "../../../src/LinearAlgebra/Models/Vector.h"
#include "../../../src/FEM/Equation/PlaneStrain.h"
#include "../../../src/FEM/Controller/ShapeFunction.h"
#include "../../../src/FEM/Controller/GaussIntegration.h"
#include "../../../src/FEM/Controller/BoundaryCondition.h"
#include "../../../src/FEM/Controller/Assembling.h"
#include "../../../src/LinearAlgebra/Solvers/CG.h"
#include "../../../src/PrePost/Export/ExportToVTK.h"
#include "../../../src/Optimize/Solver/MMA.h"
#include "../../../src/Optimize/Filter/HeavisideFilter.h"
#include "../../../src/PrePost/Mesher/SquareMesh.h"
#include "../../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    //----------Generate design region----------
	SquareMesh<double> mesh = SquareMesh<double>(60.0, 40.0, 60, 40);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 60.0) < 1.0e-5 && abs(_x(1) - 20.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -1.0;
    }

	//----------Get cg of element----------
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        cg[i] = CenterOfGravity(x, elements[i]);
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
	double E0 = 0.0001;
	double E1 = 210000.0;
	double Poisson = 0.3;
	double p = 3.0;

	double weightlimit = 0.5;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 0.5;

    MMA<double> optimizer = MMA<double>(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.2, 0.5, 0.7, 1.2, 1.0e-6);
			
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
			double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, 0.3, 1.0);
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
            double E = E1*pow(rho[i], p) + E0*(1.0 - pow(rho[i], p));
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, 0.3, 1.0);
            Vector<double> Keue = Ke*ElementVector(u, nodetoelement, elements[i]);
            Assembling(RF, Keue, nodetoglobal, nodetoelement, elements[i]);
        }

        Disassembling(r, RF, nodetoglobal);

        //--------------------Get compliance and sensitivities--------------------
        double f = scale0*std::inner_product(u.begin(), u.end(), r.begin(), 0.0);
		std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);
      
        for(int i = 0; i < elements.size(); i++){
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, 1.0, 0.3, 1.0);
			Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*(ue*(Ke*ue));
        }


        //*************************************************
        //  Filtering sensitivities
        //*************************************************
        std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);
        std::vector<double> dgds = filter.GetFilteredSensitivitis(s, dgdrho);
		
        
        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout("sample/optimize/compliance/result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(u, "u", fout, true);
        AddPointVectors(r, "r", fout, false);
		AddElementScalers(rho, "s", fout, true);
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
		optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });	
	}
	
	return 0;
}