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
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

    //----------Define parameters----------
    double a = 1.0;         //  Advection velocity
    double theta = 60.0;    //  Advection direction
    double k = 1.0e-6;      //  Diffusion coefficient

    //----------Culculate Ke Fe and Assembling----------
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
    for (auto element : elements) {
        Matrix<double> AdvectionTerm;
        Advection<double, ShapeFunction3Triangle, Gauss1Triangle>(AdvectionTerm, nodes, element, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0));
        Matrix<double> DiffusionTerm;
        Diffusion<double, ShapeFunction3Triangle, Gauss1Triangle>(DiffusionTerm, nodes, element, k);
        Matrix<double> AdvectionSUPGTerm;
        AdvectionSUPG<double, ShapeFunction3Triangle, Gauss1Triangle>(AdvectionSUPGTerm, nodes, element, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0), k);
        Matrix<double> AdvectionSCTerm;
        AdvectionShockCapturing<double, ShapeFunction3Triangle, Gauss1Triangle>(AdvectionSCTerm, nodes, element, a*cos(theta*M_PI/180.0), a*sin(theta*M_PI/180.0), k);
        Matrix<double> Ke = AdvectionTerm + DiffusionTerm + AdvectionSUPGTerm;// + AdvectionSCTerm;
        Assembling(K, Ke, element, field);
    }

    //----------Set Dirichlet Boundary Condition----------
    SetDirichlet(K, F, isufixed, ufixed, 1.0e5);

    //----------Solve System Equation----------
    CSR<double> Kmod = CSR<double>(K);
    CSR<double> M = ILU0(Kmod);
    std::vector<double> result = ILU0BiCGSTAB(Kmod, M, F, 100000, 1.0e-10);

    //----------Update T----------
    std::vector<double> T = std::vector<double>(nodes.size(), 0.0);
    FieldResultToNodeValue(result, T, field);
            
    //----------Save initial value----------
    std::ofstream fout(model_path + "result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    AddPointScalers(T, "T", fout, true);
    fout.close();
	
	return 0;
}