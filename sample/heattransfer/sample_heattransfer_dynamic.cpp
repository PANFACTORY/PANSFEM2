#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/Delaunay.h"
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
	Delaunay<double> mesher = Delaunay<double>({ { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } }, { { 3, 2, 1, 0 } }, {}, 0.1);
    std::vector<Vector<double> > x = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed0 = mesher.GenerateFixedlist({ 0 }, [](Vector<double> _x) {
        if(fabs(_x(0)) < 1.0e-5) {
            return true;
        }
        return false;
    });
	for(auto& ufixedi : ufixed0) {
		ufixedi.second = 300.0;
	}
	std::vector<std::pair<std::pair<int, int>, double> > ufixed1 = mesher.GenerateFixedlist({ 0 }, [](Vector<double> _x) {
        if(fabs(_x(0) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
	std::vector<Vector<double> > T = std::vector<Vector<double> >(x.size(), Vector<double>(1));

    std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
    SetDirichlet(T, nodetoglobal, ufixed0);
	SetDirichlet(T, nodetoglobal, ufixed1);
	int KDEGREE = Renumbering(nodetoglobal);

	double dt = 0.001;
	double theta = 0.5;

	for(int t = 0; t < 500; t++){
		std::cout << "t = " << t << std::endl;

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
		
        for (auto element : elements) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> Ke;
			HeatTransfer<double, ShapeFunction3Triangle, Gauss1Triangle>(Ke, nodetoelement, element, { 0 }, x, 1.0, 1.0);
			Matrix<double> Ce;
			HeatCapacity<double, ShapeFunction3Triangle, Gauss1Triangle>(Ce, nodetoelement, element, { 0 }, x, 1.0, 1.0, 1.0);
			Vector<double> Te = ElementVector(T, nodetoelement, element);
			Matrix<double> Ae = Ce/dt + Ke*theta;
			Vector<double> be = (Ce/dt - Ke*(1.0 - theta))*Te; 
			Assembling(K, F, T, Ae, nodetoglobal, nodetoelement, element);
            Assembling(F, be, nodetoglobal, nodetoelement, element);
		}

		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(T, result, nodetoglobal);
	
		std::ofstream fout("sample/heattransfer/result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 5), fout);
		AddPointScalers(T, "T", fout, true);
		fout.close();
	}

	return 0;
}