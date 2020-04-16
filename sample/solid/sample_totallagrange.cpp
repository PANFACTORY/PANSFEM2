#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Equation/Solid.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition2.h"
#include "../../src/FEM/Controller/Assembling2.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/solid/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
    std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "Dirichlet.csv");

	std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(3)); 
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(3, 0));
	
    SetDirichlet(u, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);

	std::ofstream fout(model_path + "result" + std::to_string(0) + ".vtk");
	MakeHeadderToVTK(fout);
	AddPointsToVTK(x, fout);
	AddElementToVTK(elements, fout);
	AddElementTypes(std::vector<int>(elements.size(), 12), fout);
	AddPointVectors(u, "u", fout, true);
	fout.close();

	for (int finc = 1, fincmax = 100; finc <= fincmax; finc++) {
		std::cout << "finc = " << finc << "\t";

		for (int k = 0; k < 100; k++) {
			LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
			std::vector<double> R = std::vector<double>(KDEGREE, 0.0);

			for (auto element : elements) {
                std::vector<std::vector<std::pair<int, int> > > nodetoelement;
				Matrix<double> Ke;
				Vector<double> Qe;
				SolidTotalLagrange<double, ShapeFunction8Cubic, Gauss8Cubic >(Ke, Qe, nodetoelement, element, { 0, 1, 2, }, x, u, 1000.0, 0.3);
				Assembling(K, R, u, Ke, Qe, nodetoglobal, nodetoelement, element);
			}

            std::vector<std::pair<std::pair<int, int>, double> > qfixed = { { { 764, 1 }, -100.0*(finc/(double)fincmax) }, { { 815, 1 }, -100.0*(finc/(double)fincmax) }, { { 1070, 1 }, -100.0*(finc/(double)fincmax) }, { { 1121, 1 }, -100.0*(finc/(double)fincmax) } };
            Assembling(R, qfixed, nodetoglobal);

			double normR = std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
            double normF = 200.0*(finc/(double)fincmax);
			if (normR/normF < 1.0e-10) {
				std::cout << "k = " << k << "\tNorm = " << normR/normF << std::endl;
				break;
			}
			
			CSR<double> Kmod = CSR<double>(K);
			std::vector<double> result = ScalingCG(Kmod, R, 100000, 1.0e-10);
			std::vector<Vector<double> > du = std::vector<Vector<double> >(x.size(), Vector<double>(3));
			Disassembling(du, result, nodetoglobal);
			
            for(int i = 0; i < x.size(); i++){
				u[i] += du[i];
			}
		}

		std::ofstream fout(model_path + "result" + std::to_string(finc) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 12), fout);
		AddPointVectors(u, "u", fout, true);
		fout.close();
	}

	return 0;
}