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
#include "../../src/FEM/Equation/Geometric.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/advection/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	std::vector<std::pair<std::pair<int, int>, double> > ufixed;
	ImportDirichletFromCSV(ufixed, model_path + "DirichletD.csv");

	std::vector<Vector<double> > T = std::vector<Vector<double> >(x.size(), Vector<double>(1));
	Vector<double> O = Vector<double>({ 0.5, 0.75 });
	for(int i = 0; i < x.size(); i++){
		double r = (x[i] - O).Norm();
		if(r <= 0.25){
			T[i](0) = 0.5*(cos(4.0*M_PI*r) + 1.0);
		} 
	}
	std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
	SetDirichlet(T, nodetoglobal, ufixed);
	int KDEGREE = Renumbering(nodetoglobal);	

	double dt = M_PI/50.0;		//	Time step
	double theta = 0.5;			//	FDM parameter for time
	double k = 0.0;				//	Diffusion coefficient

	for(int t = 0; t < 100; t++){
		std::cout << "t = " << t << std::endl;	

		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

		for (auto element : elements) {
			Vector<double> ge = CenterOfGravity(x, element);
			double ax = -(ge(1) - 0.5);
			double ay = (ge(0) - 0.5);

			Vector<double> Te = Vector<double>();
			for(auto i : element){
				Te = Te.Vstack(T[i]);
			}

			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> M;
			Mass<double, ShapeFunction3Triangle, Gauss1Triangle>(M, nodetoelement, element, { 0 }, x);
			Matrix<double> MS;
			MassSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(MS, nodetoelement, element, { 0 }, x, ax, ay, k);
			Matrix<double> A;
			Advection<double, ShapeFunction3Triangle, Gauss1Triangle>(A, nodetoelement, element, { 0 }, x, ax, ay);
			Matrix<double> D;
			Diffusion<double, ShapeFunction3Triangle, Gauss1Triangle>(D, nodetoelement, element, { 0 }, x, k);
			Matrix<double> AS;
			AdvectionSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(AS, nodetoelement, element, { 0 }, x, ax, ay, k);
			Matrix<double> Ke = (M + MS)/dt + theta*(A + D + AS);
			Vector<double> Fe = ((M + MS)/dt - (1.0 - theta)*(A + D + AS))*Te;
			Assembling(K, F, T, Ke, Fe, nodetoglobal, nodetoelement, element);
		}

		CSR<double> Kmod = CSR<double>(K);
		std::vector<double> result = BiCGSTAB(Kmod, F, 100000, 1.0e-10);
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