#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Advection.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/advection/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "ElementT.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");
	
	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "DirichletD.csv");

	//----------Initialize T----------
	std::vector<double> T = std::vector<double>(nodes.size(), 0.0);
	Vector<double> O = Vector<double>({ 0.5, 0.75 });
	for(int i = 0; i < nodes.size(); i++){
		double r = (nodes[i] - O).Norm();
		if(r <= 0.25){
			T[i] = 0.5*(cos(4.0*M_PI*r) + 1.0);
		} 
	}
		
	//----------Define time step and theta----------
	double dt = M_PI/50.0;		//	Time step
	double theta = 0.5;			//	FDM parameter for time
	double k = 0.0;				//	Diffusion coefficient

	//----------Time step loop----------
	for(int t = 0; t < 100; t++){
		std::cout << "t = " << t << std::endl;

		//----------Culculate Ke Fe and Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
		std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
		for (auto element : elements) {
			Vector<double> ge = Vector<double>(2);
			for(auto i : element){
				ge += nodes[i];
			}
			ge /= (double)element.size();

			double ax = -(ge(1) - 0.5);
			double ay = (ge(0) - 0.5);
			
			std::vector<double> Ttmp;
			for(auto i : element){
				Ttmp.push_back(T[i]);
			}
			Vector<double> Te = Vector<double>(Ttmp);

			Matrix<double> MassTerm;
			Mass<double, ShapeFunction3Triangle, Gauss1Triangle>(MassTerm, nodes, element);
			Matrix<double> MassSUPGTerm;
			MassSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(MassSUPGTerm, nodes, element, ax, ay, k);
			Matrix<double> AdvectionTerm;
			Advection<double, ShapeFunction3Triangle, Gauss1Triangle>(AdvectionTerm, nodes, element, ax, ay);
			Matrix<double> DiffusionTerm;
			Diffusion<double, ShapeFunction3Triangle, Gauss1Triangle>(DiffusionTerm, nodes, element, k);
			Matrix<double> AdvectionSUPGTerm;
			AdvectionSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(AdvectionSUPGTerm, nodes, element, ax, ay, k);

			Matrix<double> Ke = (MassTerm + MassSUPGTerm)/dt + theta*(AdvectionTerm + DiffusionTerm + AdvectionSUPGTerm);
			Vector<double> Fe = ((MassTerm + MassSUPGTerm)/dt - (1.0 - theta)*(AdvectionTerm + DiffusionTerm + AdvectionSUPGTerm))*Te;
			Assembling(K, F, Ke, Fe, element, field);
		}

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(K, F, isufixed, ufixed, 1.0e5);

		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(K);
		CSR<double> M = ILU0(Kmod);
		std::vector<double> result = ILU0BiCGSTAB(Kmod, M, F, 100000, 1.0e-10);

		//----------Update T----------
		FieldResultToNodeValue(result, T, field);
				
		//----------Save initial value----------
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 5), fout);
		AddPointScalers(T, "T", fout, true);
		fout.close();
	}

	return 0;
}