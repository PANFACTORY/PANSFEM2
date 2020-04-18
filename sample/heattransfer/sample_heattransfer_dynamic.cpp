#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/HeatTransfer.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/heattransfer/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
    std::vector<Vector<double> > T = std::vector<Vector<double> >(x.size(), Vector<double>(1));

    std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
    std::vector<std::pair<std::pair<int, int>, double> > Tfixed = { { { 0, 0 }, 100 }, { { 1, 0 }, 100 }, { { 2, 0 }, 100 }, { { 3, 0 }, 100 }, { { 4, 0 }, 100 }, { { 20, 0 }, 0 }, { { 21, 0 }, 0 }, { { 22, 0 }, 0 }, { { 23, 0 }, 0 }, { { 24, 0 }, 0 } };
    SetDirichlet(T, nodetoglobal, Tfixed);
	int KDEGREE = Renumbering(nodetoglobal);

	double dt = 0.1;
	double theta = 0.5;

	for(int t = 0; t < 100; t++){
		std::cout << "t = " << t << std::endl;

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
		
        for (auto element : elements) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> Ke;
			HeatTransfer<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, element, { 0 }, x, 1.0, 1.0);
			Matrix<double> Ce;
			HeatCapacity<double, ShapeFunction4Square, Gauss4Square>(Ce, nodetoelement, element, { 0 }, x, 1.0, 1.0, 1.0);
			Vector<double> Te = ElementVector(T, nodetoelement, element);
			Matrix<double> Ae = Ce/dt + Ke*theta;
			Vector<double> be = (Ce/dt - Ke*(1.0 - theta))*Te; 
			Assembling(K, F, T, Ae, nodetoglobal, nodetoelement, element);
            Assembling(F, be, nodetoglobal, nodetoelement, element);
		}

		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(T, result, nodetoglobal);
	
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 5), fout);
		AddPointScalers(T, "T", fout, true);
		fout.close();
	}

	return 0;
}