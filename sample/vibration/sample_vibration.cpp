#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling2.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/LinearAlgebra/Solvers/Lanczos.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"


using namespace PANSFEM2;


int main() {
	std::string model_path = "sample/bibration/";
	std::vector<Vector<double> > x;
	ImportNodesFromCSV(x, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	
    std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
	int KDEGREE = Renumbering(nodetoglobal);

    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
    LILCSR<double> M = LILCSR<double>(KDEGREE, KDEGREE);
    for (auto element : elements) {
        std::vector<std::vector<std::pair<int, int> > > nodetoelement;
		Matrix<double> Ke;
		PlaneStrainStiffness<double, ShapeFunction8Square, Gauss9Square>(Ke, nodetoelement, element, { 0, 1 }, x, 210000.0, 0.3, 1.0);
        Assembling(K, Ke, nodetoglobal, nodetoelement, element);
        Matrix<double> Me;
        PlaneStrainMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodetoelement, element, { 0, 1 }, x, 0.0000078, 1.0);       
        Assembling(M, Me, nodetoglobal, nodetoelement, element);
    }
  
    CSR<double> Kmod = CSR<double>(K);
    CSR<double> Mmod = CSR<double>(M);
    std::vector<double> eigenvalues;
	std::vector<std::vector<double> > eigenvectors;
    int m = 5;
	GeneralShiftedInvertLanczos(Kmod, Mmod, eigenvalues, eigenvectors, m, -100.0);

    for(int i = 0; i < m; i++){
        std::cout << eigenvalues[i] << "\t" << sqrt(eigenvalues[i]) << std::endl;
        
        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        Disassembling(u, eigenvectors[i], nodetoglobal);

        std::ofstream fout(model_path + "result" + std::to_string(i) + ".vtk");
        MakeHeadderToVTK(fout);
        AddPointsToVTK(x, fout);
        AddElementToVTK(elements, fout);
        AddElementTypes(std::vector<int>(elements.size(), 23), fout);
        AddPointVectors(u, "u", fout, true);
        fout.close();
    }
	
	return 0;
}