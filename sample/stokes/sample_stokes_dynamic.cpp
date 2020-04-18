#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/Stokes.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/stokes/model1/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elementsu;
	ImportElementsFromCSV(elementsu, model_path + "ElementU.csv");
	std::vector<std::vector<int> > elementsp;
	ImportElementsFromCSV(elementsp, model_path + "ElementP.csv");
	std::vector<std::vector<int> > edgesu;
	ImportElementsFromCSV(edgesu, model_path + "EdgeU.csv");
	std::vector<std::vector<int> > edgesp;
	ImportElementsFromCSV(edgesp, model_path + "EdgeP.csv");
    std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");
	
    std::vector<Vector<double> > up = std::vector<Vector<double> >(x.size(), Vector<double>(3));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));

    SetDirichlet(up, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);
    
	double dt = 0.01;   		//	Time step
	double theta = 0.5;			//	FDM parameter for time

	for(int t = 0; t < 200; t++){
		std::cout << "t = " << t << std::endl;

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
        
        for (int i = 0; i < elementsu.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelementu;
            std::vector<std::vector<std::pair<int, int> > > nodetoelementp;
            
            Vector<double> upe = Vector<double>(2*elementsu[i].size() + elementsp[i].size());
            for(int j = 0; j < elementsu[i].size(); j++){
                upe(j) = up[elementsu[i][j]](0);
                upe(j + elementsu[i].size()) = up[elementsu[i][j]](1);
            }
            for(int j = 0; j < elementsp[i].size(); j++){
                upe(j + 2*elementsu[i].size()) = up[elementsp[i][j]](2);
            }

            Matrix<double> Ke;
            Stokes<double, ShapeFunction6Triangle, ShapeFunction3Triangle, Gauss3Triangle>(Ke, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, 1.0);
            Matrix<double> Me;
            StokesMass<double, ShapeFunction6Triangle, ShapeFunction3Triangle, Gauss3Triangle>(Me, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, 1.0);
            Matrix<double> Ae = Me/dt + theta*Ke;
			Vector<double> be = (Me/dt - (1.0 - theta)*Ke)*upe;

            Assembling(K, F, up, Ae, nodetoglobal, nodetoelementu, elementsu[i], nodetoelementu, elementsu[i]);
            Assembling(K, F, up, Ae, nodetoglobal, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i]);
            Assembling(K, F, up, Ae, nodetoglobal, nodetoelementp, elementsp[i], nodetoelementu, elementsu[i]);
            Assembling(K, F, up, Ae, nodetoglobal, nodetoelementp, elementsp[i], nodetoelementp, elementsp[i]); 

            Assembling(F, be, nodetoglobal, nodetoelementu, elementsu[i]);
            Assembling(F, be, nodetoglobal, nodetoelementp, elementsp[i]);  
        }

		for(int i = 0; i < edgesu.size(); i++){
            std::vector<std::vector<std::pair<int, int> > > nodetoelementu;
            std::vector<std::vector<std::pair<int, int> > > nodetoelementp;
            Vector<double> Fe;
            StokesTraction<double, ShapeFunction3Line, ShapeFunction2Line, Gauss2Line>(Fe, nodetoelementu, edgesu[i], nodetoelementp, edgesp[i], { 0, 1, 2 }, x, [](Vector<double> _x){
                Vector<double> f = { 1.0/3.0, 0.0 };
                return f;
            });
            Assembling(F, Fe, nodetoglobal, nodetoelementu, edgesu[i]);
            Assembling(F, Fe, nodetoglobal, nodetoelementp, edgesp[i]);
        }

		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
        Disassembling(up, result, nodetoglobal);

        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
        std::vector<double> p = std::vector<double>(x.size(), 0.0);
        for(int i = 0; i < x.size(); i++){
            u[i] = up[i].Segment(0, 2);
            p[i] = up[i](2);
        }
        
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