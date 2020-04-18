#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Equation/Homogenization.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition2.h"
#include "../../src/FEM/Controller/Assembling2.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/homogenization/model2/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

    std::vector<Vector<double> > chi0 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<Vector<double> > chi1 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<Vector<double> > chi2 = std::vector<Vector<double> >(x.size(), Vector<double>(2));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
	





	
	//----------Add Dirichlet Condition----------
	std::vector<int> ismasterfixed;
	std::vector<int> isslavefixed;
	ImportPeriodicFromCSV(ismasterfixed, isslavefixed, field, model_path + "Periodic.csv");







	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F0 = std::vector<double>(KDEGREE, 0.0);
	std::vector<double> F1 = std::vector<double>(KDEGREE, 0.0);
	std::vector<double> F2 = std::vector<double>(KDEGREE, 0.0);
	
    for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, element, { 0, 1 }, x, 68.894, 0.33, 1.0);





		Assembling(K, Ke, elements[i], field);





		Matrix<double> Fes;
		HomogenizePlaneStrainBodyForce<double, ShapeFunction4Square, Gauss4Square>(Fes, nodes, elements[i], 68.894, 0.33, 1.0);
		Vector<double> Fe0 = Fes.Block(0, 0, Ke.ROW(), 1);
		Assembling(F0, Fe0, elements[i], field);
		Vector<double> Fe1 = Fes.Block(0, 1, Ke.ROW(), 1);
		Assembling(F1, Fe1, elements[i], field);
		Vector<double> Fe2 = Fes.Block(0, 2, Ke.ROW(), 1);
		Assembling(F2, Fe2, elements[i], field);
	}

	//----------Set Periodic Boundary Condition----------
	SetPeriodic(K, ismasterfixed, isslavefixed, 1.0e8);








	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);
    Disassembling(chi0, result0, nodetoglobal);
	std::vector<double> result1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
    Disassembling(chi1, result1, nodetoglobal);
	std::vector<double> result2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);
    Disassembling(chi2, result2, nodetoglobal);

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