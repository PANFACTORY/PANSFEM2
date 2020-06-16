#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/GrandStructure.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Equation/Truss.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
    //********************Set design parameters********************
    double A0 =  0.25;          //  Initial section      [mm^2]
    double E = 200000.0;        //  Young modulus        [MPa]
    double sigmat = 50.0;       //  Tensile strength     [MPa]
    double sigmac = 50.0;       //  Compressive strength [MPa]

	//********************Set model datas********************
	GrandStructure2D<double> mesh = GrandStructure2D<double>(4000.0, 2000.0, 4, 2, 1500);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5 && (abs(_x(1)) < 1.0e-5 || abs(_x(1) - 2000.0) < 1.0e-5)) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 4000.0) < 1.0e-5 && abs(_x(1) - 1000.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -5000.0;
    }
    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
    std::vector<double> A = std::vector<double>(elements.size(), A0);

	//********************Get displacement********************
    SetDirichlet(u, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

	LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			
	std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		

	for(int i = 0; i < elements.size(); i++) {
		std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		Truss2D<double>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E, A[i]);
		Assembling(K, F, u, Ke, nodetoglobal, nodetoelement, elements[i]);
	}
	
    Assembling(F, qfixed, nodetoglobal);

	CSR<double> Kmod = CSR<double>(K);
	std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
	Disassembling(u, result, nodetoglobal);

	//********************Update ********************
	for(int i = 0; i < elements.size(); i++) {
        Vector<double> li = x[elements[i][1]] - x[elements[i][0]];
        Vector<double> di = u[elements[i][1]] - u[elements[i][0]];
        double fi = E*A[i]*di*li/pow(li.Norm(), 2.0);
        if(fi > 0.0) {
            A[i] = fi/sigmat;
        } else {
            A[i] = -fi/sigmac;
        }
    }

	//********************Export datas********************
	std::ofstream fout("sample/truss/result.vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 3), fout);
	AddPointVectors(u, "u", fout, true);
	AddElementScalers(A, "A", fout, true);
	fout.close();

	return 0;
}