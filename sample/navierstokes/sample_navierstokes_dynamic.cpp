#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
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
    std::string model_path = "sample/navierstokes/model3/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elementsu;
	ImportElementsFromCSV(elementsu, model_path + "ElementU.csv");
	std::vector<std::vector<int> > elementsp;
	ImportElementsFromCSV(elementsp, model_path + "ElementP.csv");
    std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");

    std::vector<Vector<double> > up = std::vector<Vector<double> >(x.size(), Vector<double>(3));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));
	
    double rho = 1.0;
    double mu = 1.0/1000.0;
    int tmax = 200;
    double dt = 0.1;
    double theta = 0.5;
    SetDirichlet(up, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

    std::vector<Vector<double> > ubar = std::vector<Vector<double> >(x.size(), Vector<double>(2));      //  Advection velocity


    //----------Time step loop----------
    for(int t = 0; t < tmax; t++) {
        std::cout << "t=" << t << std::endl;

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//  System stiffness matrix
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//  System load vector
        
        for (int i = 0; i < elementsu.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelementu, nodetoelementp;
            Matrix<double> Ke, Me, Ce;
            NavierStokesStiffness<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Ke, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, ubar, rho, mu);
            NavierStokesConsistentMass<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Me, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, rho);
            ContinuityStiffness<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Ce, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x);
            Matrix<double> Ae = Me/dt + theta*Ke + Ce;
            Vector<double> be = (Me/dt - (1.0 - theta)*Ke)*ElementVector(up, { nodetoelementu, nodetoelementp }, { elementsu[i], elementsp[i] });
            Assembling(K, F, up, Ae, nodetoglobal, { nodetoelementu, nodetoelementp }, { elementsu[i], elementsp[i] });
            Assembling(F, be, nodetoglobal, nodetoelementu, elementsu[i]);
            Assembling(F, be, nodetoglobal, nodetoelementp, elementsp[i]);
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

        std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
        MakeHeadderToVTK(fout);
        AddPointsToVTK(x, fout);
        AddElementToVTK(elementsp, fout);
        AddElementTypes(std::vector<int>(elementsp.size(), 5), fout);
        AddPointVectors(u, "u", fout, true);
        AddPointScalers(p, "p", fout, false);
        fout.close();
    }
    
	return 0;
}