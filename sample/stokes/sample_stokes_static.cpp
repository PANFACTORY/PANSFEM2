#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Stokes.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/stokes/model1/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements for velocity u----------
	std::vector<std::vector<int> > elementsu;
	ImportElementsFromCSV(elementsu, model_path + "ElementU.csv");

    //----------Add Elements for pressure p----------
	std::vector<std::vector<int> > elementsp;
	ImportElementsFromCSV(elementsp, model_path + "ElementP.csv");

    //----------Add Edge for velocity u----------
	std::vector<std::vector<int> > edgesu;
	ImportElementsFromCSV(edgesu, model_path + "EdgeU.csv");

    //----------Add Edge for pressure p----------
	std::vector<std::vector<int> > edgesp;
	ImportElementsFromCSV(edgesp, model_path + "EdgeP.csv");

	//----------Add Field for velocity----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Initialize velocity and pressure----------
	std::vector<Vector<double> > up = std::vector<Vector<double> >(nodes.size());     //  Velocity
    
	//----------Culculate Ke Fe and Assembling----------
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);        //  System stiffness matrix
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);  //  External load vector
    for (int i = 0; i < elementsu.size(); i++) {
        Matrix<double> Ke;
        Stokes<double, ShapeFunction6Triangle, ShapeFunction3Triangle, Gauss3Triangle>(Ke, nodes, elementsu[i], elementsp[i], 1.0);
        Assembling(K, Ke, elementsu[i], field);
    }

    //----------Culculate traction term and Assembling----------
	for(int i = 0; i < edgesu.size(); i++){
		Vector<double> Fe;
		StokesTraction<double, ShapeFunction3Line, ShapeFunction2Line, Gauss2Line>(Fe, nodes, edgesu[i], edgesp[i], 1.0/3.0, 0.0);
		Assembling(F, Fe, edgesu[i], field);
	}

    //----------Set Boundary Condition----------
    SetDirichlet(K, F, isufixed, ufixed, 1.0e9);

    //----------Solve System Equation----------
    CSR<double> Kmod = CSR<double>(K);
    std::vector<double> result = CG(Kmod, F, 100000, 1.0e-20);
    FieldResultToNodeValue(result, up, field);
	std::vector<Vector<double>> u = std::vector<Vector<double>>(nodes.size());
	std::vector<double> p = std::vector<double>(nodes.size(), 0.0);
	for(int i = 0; i < nodes.size(); i++){
		u[i] = up[i].Segment(0, 2);
		if(up[i].SIZE() == 3){
			p[i] = up[i](2);
		}
	}
        
    //*************************************************		
    //  Save values
    //*************************************************
    std::ofstream fout(model_path + "result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elementsp, fout);
    AddElementTypes(std::vector<int>(elementsp.size(), 5), fout);
    AddPointVectors(u, "u", fout, true);
    AddPointScalers(p, "p", fout, false);
    fout.close();
    
	return 0;
}