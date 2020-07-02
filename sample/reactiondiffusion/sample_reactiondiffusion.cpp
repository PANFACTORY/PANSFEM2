#include <iostream>
#include <vector>
#include <random>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/ReactionDiffusion.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_real_distribution<> dist1(0.0, 1.0);

	SquareMesh<double> mesh = SquareMesh<double>(1.0, 1.0, 100, 100);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    
	std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(1));
    for(int i = 0; i < u.size(); i++) {
        u[i](0) = dist1(engine);
    }
    std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
	int KDEGREE = Renumbering(nodetoglobal);

	double dt = 0.01;
	double theta = 0.5;
    double D = 1.0;

	for(int t = 0; t < 1000; t++){
		std::cout << "t = " << t << std::endl;

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);
		
        for (int i = 0; i < elements.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> Me, Ke;
            Vector<double> Fe;
            ReactionDiffusionConsistentMass<double, ShapeFunction4Square, Gauss4Square>(Me, nodetoelement, elements[i], { 0 }, x);
            ReactionDiffusionStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0 }, x, D);
            ReactionDiffusionReaction<double, ShapeFunction4Square, Gauss4Square>(Fe, nodetoelement, elements[i], { 0 }, x, u, [&](double _u, Vector<double> _dudX) {
                return _u*(1.0 - _u);
            });
			Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
			Matrix<double> Ae = Me/dt + Ke*theta;
			Vector<double> be = (Me/dt - Ke*(1.0 - theta))*ue + Fe; 
			Assembling(K, F, u, Ae, nodetoglobal, nodetoelement, elements[i]);
            Assembling(F, be, nodetoglobal, nodetoelement, elements[i]);
		}

		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(u, result, nodetoglobal);
	
		std::ofstream fout("sample/reactiondiffusion/result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointScalers(u, "u", fout, true);
		fout.close();
	}

	return 0;
}