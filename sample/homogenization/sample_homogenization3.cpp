#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/PrePost/Mesher/SquareAnnulusMesh.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Equation/Homogenization.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/homogenization/model5/";

	SquareAnnulusMesh2<double> mesh = SquareAnnulusMesh2<double>(1.0, 1.0, 20, 20, 10, 10);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();

	/*std::vector<std::pair<int, int> > ufixed;
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
	
    for (int i = 0; i < elements.size(); i++) {
		double E = 100.0;
		double Poisson = 0.3;
		
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, Poisson, 1.0);
		Assembling(K, Ke, nodetoglobal, nodetoelement, elements[i]);
		Assembling(F0, chi0, Ke, nodetoglobal, nodetoelement, elements[i]);
		Assembling(F1, chi1, Ke, nodetoglobal, nodetoelement, elements[i]);
		Assembling(F2, chi2, Ke, nodetoglobal, nodetoelement, elements[i]);

		Matrix<double> Fes;
		HomogenizePlaneStrainBodyForce<double, ShapeFunction4Square, Gauss4Square>(Fes, nodetoelement, elements[i], { 0, 1 }, x, E, Poisson, 1.0);
		Vector<double> Fe0 = Fes.Block(0, 0, Ke.ROW(), 1);
		Assembling(F0, Fe0, nodetoglobal, nodetoelement, elements[i]);
		Vector<double> Fe1 = Fes.Block(0, 1, Ke.ROW(), 1);
		Assembling(F1, Fe1, nodetoglobal, nodetoelement, elements[i]);
		Vector<double> Fe2 = Fes.Block(0, 2, Ke.ROW(), 1);
		Assembling(F2, Fe2, nodetoglobal, nodetoelement, elements[i]);
	}

	for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		WeakSpring<double>(Ke, nodetoelement, element, { 0, 1 }, x, 1.0e-9);
		Assembling(K, Ke, nodetoglobal, nodetoelement, element);
	}

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result0 = ScalingCG(Kmod, F0, 1, 1.0e-10);
    Disassembling(chi0, result0, nodetoglobal);
	std::vector<double> result1 = ScalingCG(Kmod, F1, 1, 1.0e-10);
    Disassembling(chi1, result1, nodetoglobal);
	std::vector<double> result2 = ScalingCG(Kmod, F2, 1, 1.0e-10);
    Disassembling(chi2, result2, nodetoglobal);

	Matrix<double> CH = Matrix<double>(3, 3);
	Matrix<double> I = Matrix<double>(3, 3);
	double volume = 0.0;
	for (int i = 0; i < elements.size(); i++) {
		double E = 100.0;
		double Poisson = 0.3;
		CH += HomogenizePlaneStrainConstitutive<double, ShapeFunction4Square, Gauss4Square>(x, elements[i], chi0, chi1, chi2, E, Poisson, 1.0);
		I += HomogenizePlaneStrainCheck<double, ShapeFunction4Square, Gauss4Square>(x, elements[i], chi0, chi1, chi2, 1.0);
		volume += Area<double, ShapeFunction4Square, Gauss4Square>(x, elements[i]);
	}
	std::cout << I/volume << std::endl << CH/volume << std::endl;;
*/
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 9), fout);
	/*AddPointVectors(chi0, "chi0", fout, true);
	AddPointVectors(chi1, "chi1", fout, false);
	AddPointVectors(chi2, "chi2", fout, false);
*/	fout.close();

	return 0;
}