#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/Shell.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	SquareMesh<double> mesh = SquareMesh<double>(15.0, 5.0, 3, 1);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
	for(auto& xi : x) {
		xi = { xi(0), xi(1), 0.0 };
	}
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1, 2, 3, 4, 5 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 2 }, [](Vector<double> _x){
        if(abs(_x(0) - 15.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -10.0;
    }

	std::vector<Vector<double> > ur = std::vector<Vector<double> >(x.size(), Vector<double>(6));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(6, 0));
	
    SetDirichlet(ur, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

	for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		ShellLinearIsotropicElastic2<double, ShapeFunction4Square, Gauss1Square, Gauss2Line>(Ke, nodetoelement, element, { 0, 1, 2, 3, 4, 5 }, x, 210000.0, 0.3, 1.0);
        Assembling(K, F, ur, Ke, nodetoglobal, nodetoelement, element);
	}
    Assembling(F, qfixed, nodetoglobal);

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
    Disassembling(ur, result, nodetoglobal);

    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
    std::vector<Vector<double> > theta = std::vector<Vector<double> >(x.size());
    for(int i = 0; i < x.size(); i++){
        u[i] = ur[i].Segment(0, 3);
        theta[i] = 180.0*ur[i].Segment(3, 6)/M_PI;
    }
	
	std::ofstream fout("sample/shell/result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 9), fout);
	AddPointVectors(u, "u", fout, true);
    AddPointVectors(theta, "theta", fout, false);
	fout.close();

	return 0;
}