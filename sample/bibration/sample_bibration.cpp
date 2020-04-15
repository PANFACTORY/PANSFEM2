#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/Lanczos.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"


using namespace PANSFEM2;


int main() {
    //----------Model Path----------
	std::string model_path = "sample/bibration/";

	//----------Add Nodes----------
	std::vector<Vector<double> > X;
	ImportNodesFromCSV(X, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, X.size(), model_path + "Field.csv");
    
	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

    //----------Assembling----------
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
    LILCSR<double> M = LILCSR<double>(KDEGREE, KDEGREE);
    for (int i = 0; i < elements.size(); i++) {
        Matrix<double> Ke, Me;
        PlaneStrainStiffness<double, ShapeFunction8Square, Gauss9Square >(Ke, X, elements[i], 210000.0, 0.0, 1.0);
        PlaneStrainMass<double, ShapeFunction8Square, Gauss9Square >(Me, X, elements[i], 0.0000078, 1.0);
        Assembling(K, Ke, elements[i], field);
        Assembling(M, Me, elements[i], field);
    }
  
    //----------Set Dirichlet Boundary Condition----------
    SetDirichlet(K, isufixed, ufixed, 1.0e10);

    //----------Solve System Equation----------
    CSR<double> Kmod = CSR<double>(K);
    CSR<double> Mmod = CSR<double>(M);
    std::vector<double> eigenvalues;
	std::vector<std::vector<double> > eigenvectors;
    int m = 5;
	GeneralShiftedInvertLanczos(Kmod, Mmod, eigenvalues, eigenvectors, m, 0.0);
    for(int i = 0; i < m; i++){
        //----------Get eigen value and eigen vector----------
        std::cout << eigenvalues[i] << "\t" << sqrt(eigenvalues[i]) << std::endl;
        
        //----------Post Process----------
        std::vector<Vector<double> > u = std::vector<Vector<double> >(X.size(), Vector<double>(2));;
        FieldResultToNodeValue(eigenvectors[i], u, field);

        //----------Save file----------
        std::ofstream fout(model_path + "result" + std::to_string(i) + ".vtk");
        MakeHeadderToVTK(fout);
        AddPointsToVTK(X, fout);
        AddElementToVTK(elements, fout);
        AddElementTypes(std::vector<int>(elements.size(), 23), fout);
        AddPointVectors(u, "u", fout, true);
        fout.close();
    }
	
	return 0;
}