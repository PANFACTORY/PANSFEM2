#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/Shell.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/shell/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
    std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");
    std::vector<std::pair<std::pair<int, int>, double> > qfixed;
	ImportNeumannFromCSV(qfixed, model_path + "Neumann.csv");

    std::vector<Vector<double> > v3 = std::vector<Vector<double> >(x.size(), { 0.0, 0.0, 1.0 });    //  Director vector

	std::vector<Vector<double> > ur = std::vector<Vector<double> >(x.size(), Vector<double>(5));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(5, 0));
	
    SetDirichlet(ur, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

	for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		ShellLinearIsotropicElastic<double, ShapeFunction4Square, Gauss1Square, Gauss2Line>(Ke, nodetoelement, element, { 0, 1, 2, 3, 4 }, x, v3, 210000.0, 0.3, 1.0);
        Assembling(K, F, ur, Ke, nodetoglobal, nodetoelement, element);
	}
    Assembling(F, qfixed, nodetoglobal);

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
    Disassembling(ur, result, nodetoglobal);

    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size());
    for(int i = 0; i < x.size(); i++){
        u[i] = ur[i].Segment(0, 3);
    }
	
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 9), fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	return 0;
}