#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Mesher/SquareAnnulus.h"
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
	SquareAnnulus<double> mesh = SquareAnnulus<double>(50.0, 50.0, 150.0, 150.0, 100, 100, 100);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<int> field = mesh.GenerateFields(2);
    int KDEGREE = *(field.end() - 1);
    std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(-25.0 - 1.0e-5 < _x(0) && _x(0) < 25.0 + 1.0e-5 && -25.0 - 1.0e-5 < _x(1) && _x(1) < 25.0 + 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> ufixed = std::vector<double>(isufixed.size(), 0.0);
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
    std::vector<int> isqfixed0 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) - 75.0) < 1.0e-5 && fabs(_x(1) - 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed0 = std::vector<double>(isqfixed0.size());
    for(int i = 0; i < isqfixed0.size(); i++){
        if(i%2 == 0) {
            qfixed0[i] = -1.0/sqrt(2.0);
        } else {
            qfixed0[i] = 1.0/sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed0, qfixed0);
    std::vector<int> isqfixed1 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) + 75.0) < 1.0e-5 && fabs(_x(1) - 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed1 = std::vector<double>(isqfixed1.size());
    for(int i = 0; i < isqfixed1.size(); i++){
        if(i%2 == 0) {
            qfixed1[i] = -1.0/sqrt(2.0);
        } else {
            qfixed1[i] = -1.0/sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed1, qfixed1);
    std::vector<int> isqfixed2 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) + 75.0) < 1.0e-5 && fabs(_x(1) + 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed2 = std::vector<double>(isqfixed2.size());
    for(int i = 0; i < isqfixed2.size(); i++){
        if(i%2 == 0) {
            qfixed2[i] = 1.0/sqrt(2.0);
        } else {
            qfixed2[i] = -1.0/sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed2, qfixed2);
    std::vector<int> isqfixed3 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) - 75.0) < 1.0e-5 && fabs(_x(1) + 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed3 = std::vector<double>(isqfixed3.size());
    for(int i = 0; i < isqfixed3.size(); i++){
        if(i%2 == 0) {
            qfixed3[i] = 1.0/sqrt(2.0);
        } else {
            qfixed3[i] = 1.0/sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed3, qfixed3);




    std::vector<int> isqfixed4 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) - 75.0) < 1.0e-5 && fabs(_x(1) - 0.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed4 = std::vector<double>(isqfixed4.size());
    for(int i = 0; i < isqfixed4.size(); i++){
        if(i%2 == 0) {
            qfixed4[i] = 0.0;
        } else {
            qfixed4[i] = sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed4, qfixed4);
    std::vector<int> isqfixed5 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) - 0.0) < 1.0e-5 && fabs(_x(1) - 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed5 = std::vector<double>(isqfixed5.size());
    for(int i = 0; i < isqfixed5.size(); i++){
        if(i%2 == 0) {
            qfixed5[i] = -sqrt(2.0);
        } else {
            qfixed5[i] = 0.0;
        }
    }
    SetNeumann(F, isqfixed5, qfixed5);
    std::vector<int> isqfixed6 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) + 75.0) < 1.0e-5 && fabs(_x(1) - 0.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed6 = std::vector<double>(isqfixed6.size());
    for(int i = 0; i < isqfixed6.size(); i++){
        if(i%2 == 0) {
            qfixed6[i] = 0.0;
        } else {
            qfixed6[i] = -sqrt(2.0);
        }
    }
    SetNeumann(F, isqfixed6, qfixed6);
    std::vector<int> isqfixed7 = mesh.GenerateFixedlist(2, { 0, 1 }, [&](Vector<double> _x){
        if(fabs(_x(0) - 0.0) < 1.0e-5 && fabs(_x(1) + 75.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<double> qfixed7 = std::vector<double>(isqfixed7.size());
    for(int i = 0; i < isqfixed7.size(); i++){
        if(i%2 == 0) {
            qfixed7[i] = sqrt(2.0);
        } else {
            qfixed7[i] = 0.0;
        }
    }
    SetNeumann(F, isqfixed7, qfixed7);

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