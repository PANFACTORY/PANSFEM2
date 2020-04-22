#include <iostream>
#include <vector>


#include "../../../src/LinearAlgebra/Models/Vector.h"
#include "../../../src/LinearAlgebra/Models/Matrix.h"
#include "../../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../../src/PrePost/Import/ImportFromCSV.h"
#include "../../../src/FEM/Controller/Assembling.h"
#include "../../../src/FEM/Equation/PlaneStrain.h"
#include "../../../src/FEM/Controller/BoundaryCondition.h"
#include "../../../src/LinearAlgebra/Solvers/Lanczos.h"
#include "../../../src/PrePost/Export/ExportToVTK.h"
#include "../../../src/FEM/Controller/ShapeFunction.h"
#include "../../../src/FEM/Controller/GaussIntegration.h"
#include "../../../src/FEM/Equation/General.h"
#include "../../../src/Optimize/Solver/OC.h"
#include "../../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


int main() {
	//----------Define design parameters----------
	double E0 = 0.0001;
	double E1 = 100000.0;
	double Poisson = 0.3;
	double p = 3.0;

    double rho0 = 1.0e-12;
    double rho1 = 1.0e-8;

    int m = 10;

	double beta = 0.5;

	double scale0 = 1.0e-3;
	double scale1 = 1.0;

	double weightlimit = 0.36;

	//----------Import design region----------
	std::string model_path = "sample/optimize/vibration/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        cg[i] = CenterOfGravity(x, elements[i]);
    }

    //----------Settings for optimization----------
    std::vector<std::vector<int> > neighbors = std::vector<std::vector<int> >(elements.size(), std::vector<int>());
    std::vector<std::vector<double> > w = std::vector<std::vector<double> >(elements.size(), std::vector<double>());
    double R = 6.0;
    for(int i = 0; i < elements.size(); i++){
        for(int j = 0; j < elements.size(); j++){
            if((cg[i] - cg[j]).Norm() <= R){
                neighbors[i].push_back(j);
                w[i].push_back((R - (cg[i] - cg[j]).Norm())/R);
            }
        }
    }

    HeavisideFilter<double> filter = HeavisideFilter<double>(elements.size(), neighbors, w);
	std::vector<double> s = std::vector<double>(elements.size(), weightlimit);
	OC<double> optimizer = OC<double>(s.size(), 1.0, 0.0, 1.0e6, 1.0e-4, 0.01, std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
		
	//----------Optimize loop----------
	for(int k = 0; k < 200; k++){
		std::cout << "\nk = " << k << "\t";
        if(k%100 == 0){
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


		//**************************************************
		//	Get eigenvalue and sensitivities
		//**************************************************
        
		//----------Get eigenvalues and eigenvectors----------
		std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        
        SetDirichlet(nodetoglobal, ufixed);
        int KDEGREE = Renumbering(nodetoglobal);

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		LILCSR<double> M = LILCSR<double>(KDEGREE, KDEGREE);

		for (int i = 0; i < elements.size(); i++) {
			double Ei = E1*pow(s[i], p) + E0*(1.0 - pow(s[i], p));
            double rhoi = rho1*s[i] + rho0*(1.0 - s[i]); 
			if(i == 3974 || i == 3975){
				rhoi += 1.25e-6;
			}
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> Ke;
			PlaneStrainStiffness<double, ShapeFunction8Square, Gauss9Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, Ei, Poisson, 1.0);
			Assembling(K, Ke, nodetoglobal, nodetoelement, elements[i]);
			Matrix<double> Me;
			PlaneStrainMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodetoelement, elements[i], { 0, 1 }, x, rhoi, 1.0);       
			Assembling(M, Me, nodetoglobal, nodetoelement, elements[i]);
		}

		CSR<double> Kmod = CSR<double>(K);	
		CSR<double> Mmod = CSR<double>(M);	
        std::vector<double> eigenvalues;
		std::vector<std::vector<double> > eigenvectors;
		GeneralShiftedInvertLanczos(Kmod, Mmod, eigenvalues, eigenvectors, m, 0.0);        
        
        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        Disassembling(u, eigenvectors[0], nodetoglobal);

		//----------Get sensitivities----------
		double f = scale0*eigenvalues[0];
		std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);
		double uMu = 0.0;
      
        for(int i = 0; i < elements.size(); i++){
			double dEi = p*(E1 - E0)*pow(s[i], p - 1.0);
            double drhoi = rho1 - rho0; 
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> dKe;
            PlaneStrainStiffness<double, ShapeFunction8Square, Gauss9Square>(dKe, nodetoelement, elements[i], { 0, 1 }, x, dEi, Poisson, 1.0);
			Matrix<double> dMe;
			PlaneStrainMass<double, ShapeFunction8Square, Gauss9Square >(dMe, nodetoelement, elements[i], { 0, 1 }, x, drhoi, 1.0);       
			Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            dfdrho[i] = -scale0*ue*((dKe - eigenvalues[0]*dMe)*ue);

			double rhoi = rho1*s[i] + rho0*(1.0 - s[i]); 
			if(i == 3974 || i == 3975){
				rhoi += 1.25e-6;
			}
			Matrix<double> Me;
			PlaneStrainMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodetoelement, elements[i], { 0, 1 }, x, rhoi, 1.0);
			uMu += ue*(Me*ue);
        }

		for(int i = 0; i < elements.size(); i++) {
			dfdrho[i] /= uMu;
		}

		std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 23), fout);
		AddPointVectors(u, "u" + std::to_string(0), fout, true);
		AddElementScalers(s, "s", fout, true);
		fout.close();

		//----------Update design variables with OC method----------
		std::cout << "Objective:\t" << sqrt(f/scale0) << "\tWeight:\t" << g/scale1 << "\t";
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