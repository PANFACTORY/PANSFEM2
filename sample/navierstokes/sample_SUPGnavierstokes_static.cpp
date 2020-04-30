#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/Stokes.h"
#include "../../src/FEM/Equation/NavierStokes.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/PrePost/Mesher/Delaunay.h"


using namespace PANSFEM2;


int main() {
    Delaunay<double> mesher = Delaunay<double>({ { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } }, { { 3, 2, 1, 0 } }, {}, 0.05);
    std::vector<Vector<double> > x = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > wall0 = mesher.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x) {
        if(fabs(_x(0)) < 1.0e-5 || fabs(_x(0) - 1.0) < 1.0e-5 || fabs(_x(1)) < 1.0e-5) {
            return true;
        }
        return false;
    });
	std::vector<std::pair<std::pair<int, int>, double> > wall1 = mesher.GenerateFixedlist({ 1 }, [](Vector<double> _x) {
        if(fabs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > flow0 = mesher.GenerateFixedlist({ 0 }, [](Vector<double> _x) {
        if(fabs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& flow0i : flow0) {
		flow0i.second = 1.0;
	}

    std::vector<Vector<double> > up = std::vector<Vector<double> >(x.size(), Vector<double>(3));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));
	
    double rho = 1.0;
    double mu = 1.0/100.0;
    int kmax = 10;
    
    //----------Incremental step loop----------
    for(int k = 0; k < kmax; k++) {
        std::cout << "k=" << k;

        SetDirichlet(up, nodetoglobal, wall0);
        SetDirichlet(up, nodetoglobal, wall1);
        std::vector<std::pair<std::pair<int, int>, double> > flowk = flow0;
        for(auto& flowki : flowk) {
            flowki.second *= (k + 1)/(double)kmax;
        }
        SetDirichlet(up, nodetoglobal, flowk);
        int KDEGREE = Renumbering(nodetoglobal);

        //----------Newton-Raphson loop----------
        for(int l = 0; l < 100; l++) {
            LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
            std::vector<double> R = std::vector<double>(KDEGREE, 0.0);		//Residual load vector
            std::vector<Vector<double> > dup = std::vector<Vector<double> >(x.size(), Vector<double>(3));

            for (auto element : elements) {
                std::vector<std::vector<std::pair<int, int> > > nodetoelementu, nodetoelementp;
                Matrix<double> Ke;
                Vector<double> Qe;
                NavierStokesPSPGSUPGTangent<double, ShapeFunction3Triangle, ShapeFunction3Triangle, Gauss1Triangle>(Ke, Qe, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, up, rho, mu);
                Assembling(K, R, dup, Ke, nodetoglobal, { nodetoelementu, nodetoelementp }, { element, element });
                Assembling(R, Qe, nodetoglobal, nodetoelementu, element);
                Assembling(R, Qe, nodetoglobal, nodetoelementp, element);  
            }
  
            double normR = std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
            std::cout << "\t\tl = " << l << "\tR Norm = " << normR << std::endl;
            if (normR < 1.0e-10) {
                std::cout << "\tConvergence at l = " << l << "\tR Norm = " << normR << std::endl;
                break;
            }
            
            CSR<double> Kmod = CSR<double>(K);
            std::vector<double> result = BiCGSTAB2(Kmod, R, 100000, 1.0e-10);
            Disassembling(dup, result, nodetoglobal);
            for(int i = 0; i < x.size(); i++){
                up[i] += dup[i];
            }
        }
    }

    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
    std::vector<double> p = std::vector<double>(x.size(), 0.0);
    for(int i = 0; i < x.size(); i++){
        u[i] = up[i].Segment(0, 2);
        p[i] = up[i](2);
    }

    std::ofstream fout("sample/navierstokes/SUPGresult.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(x, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    AddPointVectors(u, "u", fout, true);
    AddPointScalers(p, "p", fout, false);
    fout.close();
    
	return 0;
}