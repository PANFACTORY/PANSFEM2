#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/Assembling2.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::vector<Vector<double> > x = { { 0, 0 }, { 1, 0 }, { 2, 0 }, { 2, 1 }, { 1, 1 }, { 0, 1 } };
	std::vector<std::vector<int> > elements = { { 0, 1, 4 }, { 1, 2, 3 }, { 3, 4, 1 }, { 4, 5, 0 } };
	std::vector<std::vector<int> > edges = { { 3, 4 }, { 4, 5 } };

	std::vector<std::vector<int> > nodetoglobal = { { -1, -1 }, { 1, 2 }, { 0, 3 }, { 4, 5 }, { 6, 7 }, { -1, 8 } };
	int KDEGREE = Renumbering(nodetoglobal);

	std::vector<Vector<double> > u = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
	std::vector<Vector<double> > f = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, -100 }, { 0, 0 }, { 0, 0 } };
	
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		

	ConvertNodeToGlobal(F, f, nodetoglobal);

	for(auto element : elements) {
		std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		PlaneStrainStiffness<double, ShapeFunction3Triangle, Gauss1Triangle>(Ke, nodetoelement, element, { 0, 1 }, x, 210000.0, 0.3, 1.0);
		Assembling(K, F, u, Ke, nodetoglobal, nodetoelement, element);
		Vector<double> Fe;
		PlaneStrainBodyForce<double, ShapeFunction3Triangle, Gauss1Triangle>(Fe, nodetoelement, element, { 0, 1 }, x, [](Vector<double> _x){
			Vector<double> f = { 0.0, -300.0 };
			return f;
		}, 1.0);
		Assembling(F, Fe, nodetoglobal, nodetoelement, element);
	}

	for(auto edge : edges) {
		std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Vector<double> Fe;
		PlaneStrainSurfaceForce<double, ShapeFunction2Line, Gauss1Line>(Fe, nodetoelement, edge, { 0, 1 }, x, [](Vector<double> _x){
			Vector<double> f = { 0.0, -200.0 };
			return f;
		}, 1.0);
		Assembling(F, Fe, nodetoglobal, nodetoelement, edge);
	}

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
	Disassembling(u, result, nodetoglobal);

	std::ofstream fout("sample/planestrain/result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 5), fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	return 0;
}