#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/NavierStokes.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    SquareMesh<double> mesher = SquareMesh<double>(1, 1, 30, 30);
    std::vector<Vector<double> > x = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed0 = mesher.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x) {
        if(fabs(_x(0) - 0.0) < 1.0e-5 || fabs(_x(1) - 0.0) < 1.0e-5 || fabs(_x(0) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > ufixed1 = mesher.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x) {
        if(fabs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& ufixed1i : ufixed1) {
        if(ufixed1i.first.second == 0) {
            ufixed1i.second = 1.0;
        }
    }

    std::vector<Vector<double> > up = std::vector<Vector<double> >(x.size(), Vector<double>(3));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));
	
    double rho = 1.0;
    double mu = 1.0/1000.0;
    int tmax = 10000;
    double dt = 0.01;
    double theta = 0.5;
    SetDirichlet(up, nodetoglobal, ufixed0);
    SetDirichlet(up, nodetoglobal, ufixed1);
	int KDEGREE = Renumbering(nodetoglobal);

    std::vector<Vector<double> > ubar = std::vector<Vector<double> >(x.size(), Vector<double>(2));      //  Advection velocity


    //----------Time step loop----------
    for(int t = 0; t <= tmax; t++) {
        std::cout << "t=" << t << std::endl;

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//  System stiffness matrix
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//  System load vector
        
        for (auto element : elements) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelementu, nodetoelementp;
            Matrix<double> Ke, Me, Ce, Ks, Ms, Cs;
            NavierStokesStiffness<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, ubar, rho, mu);
            NavierStokesConsistentMass<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Me, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, rho);
            ContinuityStiffness<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Ce, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x);
            NavierStokesSUPGPSPGStiffness<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Ks, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, ubar, rho, mu, dt);
            NavierStokesSUPGPSPGConsistentMass<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Ms, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, ubar, rho, mu, dt);
            ContinuitySUPGPSPGStiffness<double, ShapeFunction4Square, ShapeFunction4Square, Gauss4Square>(Cs, nodetoelementu, element, nodetoelementp, element, { 0, 1, 2 }, x, ubar, rho, mu, dt);
            Matrix<double> Ae = (Me + Ms)/dt + theta*(Ke + Ks) + (Ce + Cs);
            Vector<double> be = ((Me + Ms)/dt - (1.0 - theta)*(Ke + Ks))*ElementVector(up, { nodetoelementu, nodetoelementp }, { element, element });
            Assembling(K, F, up, Ae, nodetoglobal, { nodetoelementu, nodetoelementp }, { element, element });
            Assembling(F, be, nodetoglobal, nodetoelementu, element);
            Assembling(F, be, nodetoglobal, nodetoelementp, element);
        }

        CSR<double> Kmod = CSR<double>(K);
        std::vector<double> result = BiCGSTAB2(Kmod, F, 100000, 1.0e-10);
        Disassembling(up, result, nodetoglobal);

        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
        std::vector<double> p = std::vector<double>(x.size(), 0.0);
        for(int i = 0; i < x.size(); i++){
            u[i] = up[i].Segment(0, 2);
            p[i] = up[i](2);
        }

        ubar = u;

        if(t%100 == 0) {
            std::ofstream fout("sample/navierstokes/result" + std::to_string(t/100) + ".vtk");
            MakeHeadderToVTK(fout);
            AddPointsToVTK(x, fout);
            AddElementToVTK(elements, fout);
            AddElementTypes(std::vector<int>(elements.size(), 5), fout);
            AddPointVectors(u, "u", fout, true);
            AddPointScalers(p, "p", fout, false);
            fout.close();
        }
    }
    
	return 0;
}