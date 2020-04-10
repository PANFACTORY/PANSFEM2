#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Mesher/Annulus.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Equation/Geometric.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/Optimize/Solver/OC.h"
#include "../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


int main() {
    //----------Generate design region----------
	Annulus<double> mesh = Annulus<double>(50.0, 150.0, 100, 360);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<int> field = mesh.GenerateFields(2);
    int KDEGREE = *(field.end() - 1);
    std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(_x.Norm() - 50.0 < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> ufixed = std::vector<double>(isufixed.size(), 0.0);
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
    int nt = 12;
    for(int i = 0; i < nt; i++) {
        double angle0 = 2.0*M_PI*(i + 1)/(double)nt;
        std::vector<int> isqfixed = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
            double angle = atan2(_x(1), _x(0)) + M_PI;
            if(fabs(angle - angle0) < 1.0e-5 && fabs(_x.Norm() - 150.0) < 1.0e-5) {
                return true;
            }
            return false;
        });
        std::vector<double> qfixed = std::vector<double>(isqfixed.size());
        for(int i = 0; i < isqfixed.size(); i++){
            if(i%2 == 0) {
                qfixed[i] = sin(angle0);
            } else {
                qfixed[i] = -cos(angle0);
            }
        }
        SetNeumann(F, isqfixed, qfixed);
    }

	//----------Get cg of element----------
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        cg[i] = CenterOfGravity(nodes, elements[i]);
    }

    //----------Get neighbor elements list and weight----------
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

    //----------Get element volume----------
    std::vector<double> areas = std::vector<double>(elements.size());
    double wholearea = 0.0;
    for(int i = 0; i < elements.size(); i++){
        areas[i] = Area<double, ShapeFunction4Square, Gauss4Square >(nodes, elements[i]);
        wholearea += areas[i];
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

	double weightlimit = 0.4;
	double scale0 = 1.0e5;
	double scale1 = 1.0;

    double beta = 0.5;

    OC<double> optimizer = OC<double>(s.size(), 0.5, 0.0, 1.0e4, 1.0e-3, 0.15, std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
			
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
            g += scale1*rho[i]*areas[i]/(weightlimit*wholearea);
            dgdrho[i] = scale1*areas[i]/(weightlimit*wholearea); 
        }
        g -= 1.0*scale1;

        
        //*************************************************
        //  Get compliance value and sensitivities
        //*************************************************
        double f = 0.0;													    //Function value of compliance
		std::vector<double> dfdrho = std::vector<double>(s.size(), 0.0);    //Sensitivities of compliance

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
            dfdrho[i] = -scale0*p*(- E0 + E1)*pow(rho[i], p - 1.0)*(de*(Ke*de));
        }


        //*************************************************
        //  Filtering sensitivities
        //*************************************************
        std::vector<double> dfds = filter.GetFilteredSensitivitis(s, dfdrho);
        std::vector<double> dgds = filter.GetFilteredSensitivitis(s, dgdrho);
		
        
        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout("sample/optimize/result" + std::to_string(k) + ".vtk");
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
        std::cout << "Objective:\t" << f/scale0 << "\tWeight:\t" << g/scale1 << "\t";
		if(k > 200 && optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with MMA----------
		optimizer.UpdateVariables(s, f, dfds, g, dgds, 
            [&](std::vector<double> _xkp1) {
                double g = 0.0;
                std::vector<double> rho = filter.GetFilteredVariables(_xkp1);
                for(int i = 0; i < elements.size(); i++){
                    g += scale1*rho[i]*areas[i]/(weightlimit*wholearea);
                }
                return g - 1.0*scale1; 
            }
        );	
	}
	
	return 0;
}