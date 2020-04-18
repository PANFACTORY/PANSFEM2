#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Equation/Homogenization.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/homogenization/model2/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
    std::vector<std::pair<int, int> > ufixed;
	ImportPeriodicFromCSV(ufixed, model_path + "Periodic.csv");

    std::vector<Vector<double> > chi0 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<Vector<double> > chi1 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<Vector<double> > chi2 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
	int KDEGREE = SetPeriodic(nodetoglobal, ufixed);

	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F0 = std::vector<double>(KDEGREE, 0.0);
	std::vector<double> F1 = std::vector<double>(KDEGREE, 0.0);
	std::vector<double> F2 = std::vector<double>(KDEGREE, 0.0);
	
    for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, element, { 0, 1 }, x, 68.894, 0.33, 1.0);
		Assembling(K, Ke, nodetoglobal, nodetoelement, element);
		Assembling(F0, chi0, Ke, nodetoglobal, nodetoelement, element);
		Assembling(F1, chi1, Ke, nodetoglobal, nodetoelement, element);
		Assembling(F2, chi2, Ke, nodetoglobal, nodetoelement, element);

		Matrix<double> Fes;
		HomogenizePlaneStrainBodyForce<double, ShapeFunction4Square, Gauss4Square>(Fes, nodetoelement, element, { 0, 1 }, x, 68.894, 0.33, 1.0);
		Vector<double> Fe0 = Fes.Block(0, 0, Ke.ROW(), 1);
		Assembling(F0, Fe0, nodetoglobal, nodetoelement, element);
		Vector<double> Fe1 = Fes.Block(0, 1, Ke.ROW(), 1);
		Assembling(F1, Fe1, nodetoglobal, nodetoelement, element);
		Vector<double> Fe2 = Fes.Block(0, 2, Ke.ROW(), 1);
		Assembling(F2, Fe2, nodetoglobal, nodetoelement, element);
	}

	for(int i = 0; i < KDEGREE; i++) {
		K.set(i, i, K.get(i, i) + 1.0e-6);
	}

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);
    Disassembling(chi0, result0, nodetoglobal);
	std::vector<double> result1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
    Disassembling(chi1, result1, nodetoglobal);
	std::vector<double> result2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);
    Disassembling(chi2, result2, nodetoglobal);

	Matrix<double> CH = Matrix<double>(3, 3);
	for (auto element : elements) {
		CH += HomogenizePlaneStrainConstitutive<double, ShapeFunction4Square, Gauss4Square>(x, element, chi0, chi1, chi2, 68.894, 0.33, 1.0);
	}
	std::cout << CH;

	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 9), fout);
	AddPointVectors(chi0, "chi0", fout, true);
	AddPointVectors(chi1, "chi1", fout, false);
	AddPointVectors(chi2, "chi2", fout, false);
	fout.close();

	return 0;
}