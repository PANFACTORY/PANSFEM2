#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/Advection.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/advection/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
    std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");

    double a = 1.0;         //  Advection velocity
    double theta = 60.0;    //  Advection direction
    double k = 1.0e-6;      //  Diffusion coefficient

    std::vector<Vector<double> > T = std::vector<Vector<double> >(x.size(), Vector<double>(1));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
	
    SetDirichlet(T, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

    for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
        Matrix<double> A;
        Advection<double, ShapeFunction3Triangle, Gauss1Triangle>(A, nodetoelement, element, { 0 }, x, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0));
        Matrix<double> B;
        Diffusion<double, ShapeFunction3Triangle, Gauss1Triangle>(B, nodetoelement, element, { 0 }, x, k);
        Matrix<double> C;
        AdvectionSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(C, nodetoelement, element, { 0 }, x, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0), k);
        Matrix<double> D;
        AdvectionShockCapturing<double, ShapeFunction3Triangle, Gauss1Triangle>(D, nodetoelement, element, { 0 }, x, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0), k);
        Matrix<double> Ke = A + B + C;// + AdvectionSCTerm;
        Assembling(K, F, T, Ke, nodetoglobal, nodetoelement, element);
    }

    CSR<double> Kmod = CSR<double>(K);
    std::vector<double> result = BiCGSTAB(Kmod, F, 100000, 1.0e-10);
    Disassembling(T, result, nodetoglobal);
            
    std::ofstream fout(model_path + "result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(x, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    AddPointScalers(T, "T", fout, true);
    fout.close();
	
	return 0;
}