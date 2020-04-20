#include <iostream>
#include <vector>


#include "../../../src/LinearAlgebra/Models/Vector.h"
#include "../../../src/FEM/Equation/PlaneStrain.h"
#include "../../../src/FEM/Controller/ShapeFunction.h"
#include "../../../src/FEM/Controller/GaussIntegration.h"
#include "../../../src/FEM/Controller/BoundaryCondition.h"
#include "../../../src/FEM/Controller/Assembling.h"
#include "../../../src/LinearAlgebra/Solvers/CG.h"
#include "../../../src/PrePost/Mesher/SquareMesh.h"
#include "../../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    //----------Define design parameters----------
	double E1 = 210000.0;
	double Poisson = 0.3;

    //----------Generate design region----------
	SquareMesh<double> mesh = SquareMesh<double>(100.0, 100.0, 100, 100);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(1)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 20.0) < 1.0e-5 && abs(_x(1) - 100.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -1.0;
    }
	
    //--------------------Get displacement--------------------
    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
    
    SetDirichlet(u, nodetoglobal, ufixed);
    int KDEGREE = Renumbering(nodetoglobal);

    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

    for (int i = 0; i < elements.size(); i++) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
        Matrix<double> Ke;
        PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E1, Poisson, 1.0);
        Assembling(K, Ke, nodetoglobal, nodetoelement, elements[i]);
        Assembling(F, u, Ke, nodetoglobal, nodetoelement, elements[i]);
    }
    Assembling(F, qfixed, nodetoglobal);

    CSR<double> Kmod = CSR<double>(K);	
    std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
    Disassembling(u, result, nodetoglobal);

    //--------------------Get reaction force--------------------
    RemoveBoundaryConditions(nodetoglobal);
    KDEGREE = Renumbering(nodetoglobal);
    F = std::vector<double>(KDEGREE, 0.0);
    std::vector<Vector<double> > r = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    
    for (int i = 0; i < elements.size(); i++) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
        Matrix<double> Ke;
        PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, E1, Poisson, 1.0);
        Vector<double> Keue = Ke*ElementVector(u, nodetoelement, elements[i]);
        Assembling(F, Keue, nodetoglobal, nodetoelement, elements[i]);
    }

    Disassembling(r, F, nodetoglobal);

    //--------------------Get compliance and sensitivities--------------------
    std::cout << std::inner_product(r.begin(), r.end(), u.begin(), 0.0) << std::endl;;
	
	return 0;
}