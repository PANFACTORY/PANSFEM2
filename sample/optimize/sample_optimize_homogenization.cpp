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
#include "../../src/PrePost/Mesher/SquareAnnulusMesh.h"
#include "../../src/FEM/Equation/General.h"
#include "../../src/FEM/Equation/Homogenization.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"


using namespace PANSFEM2;


int main() {
    //*****************************************************
    //  Define design parameters
    //*****************************************************
    double E0 = 100.0;      //	Young modulus of an unit cell
	double Poisson0 = 0.3;  //	Poisson ratio of an unit cell
	double Volume0 = 1.0;   //	Volume of an unit cell

    std::string model_path = "sample/optimize/";

    std::vector<double> as = { 1.0e-3, 0.2, 0.4, 0.6, 0.8, 0.999 };

    double weightlimit = 0.5;
	

    //*****************************************************
    //  Numerical Material Experiment
    //*****************************************************
    std::vector<Matrix<double> > CH = std::vector<Matrix<double> > (as.size()); 
    for(int i = 0; i < as.size(); i++) {
        //----------Make model of an unit cell----------
        SquareAnnulusMesh<double> mesh = SquareAnnulusMesh<double>(1.0, 1.0, 1.0 - as[i], 1.0 - as[i], 10, 10, 10);
        std::vector<Vector<double> > x = mesh.GenerateNodes();
        std::vector<std::vector<int> > elements = mesh.GenerateElements();
        std::vector<std::pair<int, int> > ufixed;
        ImportPeriodicFromCSV(ufixed, model_path + "Periodic.csv");

        std::vector<Vector<double> > chi0 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > chi1 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<Vector<double> > chi2 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        int KDEGREE = SetPeriodic(nodetoglobal, ufixed);


        //----------Get characteristics displacement----------
        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F0 = std::vector<double>(KDEGREE, 0.0);
        std::vector<double> F1 = std::vector<double>(KDEGREE, 0.0);
        std::vector<double> F2 = std::vector<double>(KDEGREE, 0.0);
        
        for (auto element : elements) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, element, { 0, 1 }, x, E0, Poisson0, 1.0);
            Assembling(K, Ke, nodetoglobal, nodetoelement, element);
            Assembling(F0, chi0, Ke, nodetoglobal, nodetoelement, element);
            Assembling(F1, chi1, Ke, nodetoglobal, nodetoelement, element);
            Assembling(F2, chi2, Ke, nodetoglobal, nodetoelement, element);

            Matrix<double> Fes;
            HomogenizePlaneStrainBodyForce<double, ShapeFunction4Square, Gauss4Square>(Fes, nodetoelement, element, { 0, 1 }, x, E0, Poisson0, 1.0);
            Vector<double> Fe0 = Fes.Block(0, 0, Ke.ROW(), 1);
            Assembling(F0, Fe0, nodetoglobal, nodetoelement, element);
            Vector<double> Fe1 = Fes.Block(0, 1, Ke.ROW(), 1);
            Assembling(F1, Fe1, nodetoglobal, nodetoelement, element);
            Vector<double> Fe2 = Fes.Block(0, 2, Ke.ROW(), 1);
            Assembling(F2, Fe2, nodetoglobal, nodetoelement, element);
        }

        for (auto element : elements) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            WeakSpring<double>(Ke, nodetoelement, element, { 0, 1 }, x, 1.0e-9);
            Assembling(K, Ke, nodetoglobal, nodetoelement, element);
        }

        CSR<double> Kmod = CSR<double>(K);
        std::vector<double> result0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);
        Disassembling(chi0, result0, nodetoglobal);
        std::vector<double> result1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
        Disassembling(chi1, result1, nodetoglobal);
        std::vector<double> result2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);
        Disassembling(chi2, result2, nodetoglobal);


        //----------Get homogenized value----------
        CH[i] = Matrix<double>(3, 3);
        for (auto element : elements) {
            CH[i] += HomogenizePlaneStrainConstitutive<double, ShapeFunction4Square, Gauss4Square>(x, element, chi0, chi1, chi2, E0, Poisson0, 1.0);
        }
        CH[i] /= Volume0;

        //----------Export result of Numerical Material Experiment----------
        std::ofstream fout(model_path + "result_microscopic_a" + std::to_string(i) + ".vtk");
        MakeHeadderToVTK(fout);
        AddPointsToVTK(x, fout);
        AddElementToVTK(elements, fout);
        AddElementTypes(std::vector<int>(elements.size(), 9), fout);
        AddPointVectors(chi0, "chi0", fout, true);
        AddPointVectors(chi1, "chi1", fout, false);
        AddPointVectors(chi2, "chi2", fout, false);
        fout.close();
    }

    
    //*****************************************************
    //  Define Macroscopic problem
    //*****************************************************
    
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
        if(abs(_x(0) - 60.0) < 1.0e-5 && abs(_x(1) - 20.0) < 2.0 + 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -1.0;
    }


	//----------Initialize design variables and solver----------
	std::vector<double> a = std::vector<double>(elements.size(), 0.5);
    OC<double> optimizer = OC<double>(a.size(), 1.0, 0.0, 1.0e4, 1.0e-3, 0.05, std::vector<double>(a.size(), 0.01), std::vector<double>(a.size(), 1.0));
			

	//----------Optimize loop----------
	for(int k = 0; k < 400; k++){
		std::cout << "\nk = " << k << "\t";
        

        //*************************************************
        //  Get weight value and sensitivities
        //*************************************************
        double g = 0.0;														//Function values of weight
		std::vector<double> dgda = std::vector<double>(a.size(), 0.0);      //Sensitivities of weight
        for(int i = 0; i < elements.size(); i++){
            g += (1.0 - pow(1.0 - a[i], 2.0))/(weightlimit*elements.size());
            dgda[i] = 2.0*(1.0 - a[i])/(weightlimit*elements.size()); 
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
			Matrix<double> CHi = Matrix<double>(3, 3);  
            for(int n = 0; n < as.size(); n++) {
                Matrix<double> P = CH[n];
                for(int m = 0; m < as.size(); m++) {
                    if(m != n) {
                        P *= (a[i] - as[m])/(as[n] - as[m]);
                    }
                }
                CHi += P;
            }
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, CHi, 1.0);
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
            Matrix<double> CHi = Matrix<double>(3, 3);
            for(int n = 0; n < as.size(); n++) {
                Matrix<double> P = CH[n];
                for(int m = 0; m < as.size(); m++) {
                    if(m != n) {
                        P *= (a[i] - as[m])/(as[n] - as[m]);
                    }
                }
                CHi += P;
            }      
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, CHi, 1.0);
            Vector<double> Keue = Ke*ElementVector(u, nodetoelement, elements[i]);
            Assembling(RF, Keue, nodetoglobal, nodetoelement, elements[i]);
        }

        Disassembling(r, RF, nodetoglobal);

        //--------------------Get compliance and sensitivities--------------------
        double f = std::inner_product(u.begin(), u.end(), r.begin(), 0.0);
		std::vector<double> dfda = std::vector<double>(a.size(), 0.0);
      
        for(int i = 0; i < elements.size(); i++){
            Matrix<double> dCHida = Matrix<double>(3, 3);
            for(int n = 0; n < as.size(); n++) {
                double q = 0.0;
                for(int m = 0; m < as.size(); m++) {
                    if(m != n ) {
                        double p = 1.0;
                        for(int l = 0; l < as.size(); l++) {
                            if(l != m && l != n) {
                                p *= (a[i] - as[l])/(as[n] - as[l]);
                            } else if(l != n) {
                                p *= 1.0/(as[n] - as[l]); 
                            }
                        }
                        q += p;
                    }
                }
                dCHida += q*CH[n];
            }
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, dCHida, 1.0);
			Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            dfda[i] = -ue*(Ke*ue);
        }
		
         
        //*************************************************
        //  Export result of Macroscopic optimization
        //*************************************************
		std::ofstream fout(model_path + "result_macroscopic_k" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(u, "u", fout, true);
        AddPointVectors(r, "r", fout, false);
		AddElementScalers(a, "a", fout, true);
		fout.close();
       

        //*************************************************
        //  Update design variable a
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << f << "\tWeight:\t" << g << "\t";
		if(optimizer.IsConvergence(f)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with OC----------
		optimizer.UpdateVariables(a, f, dfda, g, dgda, 
            [&](std::vector<double> _a) {
                double g = 0.0;
                for(int i = 0; i < _a.size(); i++){
                    g += (1.0 - pow(1.0 - _a[i], 2.0))/(weightlimit*elements.size());
                }
                return g - 1.0; 
            }
        );	
	}
	
	return 0;
}