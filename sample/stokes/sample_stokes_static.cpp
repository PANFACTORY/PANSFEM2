#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/Stokes.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	SquareMesh<double> mesh = SquareMesh<double>(1.0, 1.0, 20, 20);
    std::vector<Vector<double> > x = mesh.GenerateNodes2();
    std::vector<std::vector<int> > elementsu = mesh.GenerateElements2();
	std::vector<std::vector<int> > elementsp = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed0 = mesh.GenerateFixedlist2({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5 || abs(_x(1)) < 1.0e-5 || abs(_x(0) - 1.0) < 1.0e-5 || abs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
	std::vector<std::pair<std::pair<int, int>, double> > ufixed1 = mesh.GenerateFixedlist2({ 0 }, [](Vector<double> _x){
        if(abs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
	for(auto& ufixed1i : ufixed1) {
		ufixed1i.second = 1.0;
	}
	
    std::vector<Vector<double> > up = std::vector<Vector<double> >(x.size(), Vector<double>(3));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));
	
	SetDirichlet(up, nodetoglobal, ufixed0);
	SetDirichlet(up, nodetoglobal, ufixed1);
	int KDEGREE = Renumbering(nodetoglobal);
    
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE); 
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

    for (int i = 0; i < elementsu.size(); i++) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelementu;
        std::vector<std::vector<std::pair<int, int> > > nodetoelementp;
        Matrix<double> Ke;
        StokesStiffness<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Ke, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, 1.0);
        Assembling(K, F, up, Ke, nodetoglobal, { nodetoelementu, nodetoelementp }, { elementsu[i], elementsp[i] });
		Vector<double> Fe;
		StokesBodyForce<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Fe, nodetoelementu, elementsu[i], nodetoelementp, elementsp[i], { 0, 1, 2 }, x, [](Vector<double> _x){
			Vector<double> f = { 0.0, 0.0};
			if(_x(0) > 0.5 && _x(1) < 0.5) {
				f(1) = -100.0;
			}
			return f;
		});
		Assembling(F, Fe, nodetoglobal, nodetoelementu, elementsu[i]);
    }

	/*for(int i = 0; i < edgesu.size(); i++){
        std::vector<std::vector<std::pair<int, int> > > nodetoelementu;
        std::vector<std::vector<std::pair<int, int> > > nodetoelementp;
		Vector<double> Fe;
		StokesSurfaceForce<double, ShapeFunction3Line, ShapeFunction2Line, Gauss2Line>(Fe, nodetoelementu, edgesu[i], nodetoelementp, edgesp[i], { 0, 1, 2 }, x, [](Vector<double> _x){
			Vector<double> f = { 1.0/3.0, 0.0 };
			return f;
		});
		Assembling(F, Fe, nodetoglobal, nodetoelementu, edgesu[i]);
        Assembling(F, Fe, nodetoglobal, nodetoelementp, edgesp[i]);
	}*/

    CSR<double> Kmod = CSR<double>(K);
    std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
    Disassembling(up, result, nodetoglobal);

	std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
	std::vector<double> p = std::vector<double>(x.size(), 0.0);
	for(int i = 0; i < x.size(); i++){
		u[i] = up[i].Segment(0, 2);
		p[i] = up[i](2);
	}
        
    std::ofstream fout("sample/stokes/result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(x, fout);
    AddElementToVTK(elementsp, fout);
    AddElementTypes(std::vector<int>(elementsp.size(), 9), fout);
    AddPointVectors(u, "u", fout, true);
    AddPointScalers(p, "p", fout, false);
    fout.close();
    
	return 0;
}