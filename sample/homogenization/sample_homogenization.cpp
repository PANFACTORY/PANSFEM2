#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/homogenization/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");
	
	//----------Add Dirichlet Condition----------
	std::vector<int> ismasterfixed;
	std::vector<int> isslavefixed;
	ImportPeriodicFromCSV(ismasterfixed, isslavefixed, field, model_path + "Periodic.csv");

	//----------Culculate Ke, body force and Assembling----------
	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
	std::vector<double> F0 = std::vector<double>(KDEGREE, 0.0);		//External load vector
	std::vector<double> F1 = std::vector<double>(KDEGREE, 0.0);		//External load vector
	std::vector<double> F2 = std::vector<double>(KDEGREE, 0.0);		//External load vector
	for (auto element : elements) {
		Matrix<double> Ke;
		PlaneStrain<double, ShapeFunction4Square, Gauss4Square>(Ke, nodes, element, 210000.0, 0.3, 1.0);
		Assembling(K, Ke, element, field);
		Matrix<double> Fes;
		PlaneHomogenizeForce<double, ShapeFunction4Square, Gauss4Square>(Fes, nodes, element, 210000.0, 0.3, 1.0);
		Vector<double> Fe0 = Fes.Block(0, 0, Ke.ROW(), 1);
		Assembling(F0, Fe0, element, field);
		Vector<double> Fe1 = Fes.Block(0, 1, Ke.ROW(), 1);
		Assembling(F1, Fe1, element, field);
		Vector<double> Fe2 = Fes.Block(0, 2, Ke.ROW(), 1);
		Assembling(F2, Fe2, element, field);
	}

	//----------Set Periodic Boundary Condition----------
	SetPeriodic(K, ismasterfixed, isslavefixed, 1e-7);

	//----------Solve System Equation----------
	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result0 = ScalingCG(Kmod, F0, 100000, 1.0e-10);
	std::vector<double> result1 = ScalingCG(Kmod, F1, 100000, 1.0e-10);
	std::vector<double> result2 = ScalingCG(Kmod, F2, 100000, 1.0e-10);

	//----------Update displacement u----------
	std::vector<Vector<double> > chi0 = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
	FieldResultToNodeValue(result0, chi0, field);
	std::vector<Vector<double> > chi1 = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
	FieldResultToNodeValue(result1, chi1, field);
	std::vector<Vector<double> > chi2 = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));
	FieldResultToNodeValue(result2, chi2, field);
			
	//----------Save initial value----------
	std::ofstream fout(model_path + "result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(nodes, fout);
	AddElementToVTK(elements, fout);
	std::vector<int> et = std::vector<int>(elements.size(), 9);
	AddElementTypes(et, fout);
	AddPointVectors(chi0, "chi0", fout, true);
	AddPointVectors(chi1, "chi1", fout, false);
	AddPointVectors(chi2, "chi2", fout, false);
	fout.close();

	return 0;
}