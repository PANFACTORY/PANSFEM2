#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Stokes.h"
#include "../../src/FEM/Equation/NavierStokes.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"


using namespace PANSFEM2;


int main() {
    //*****************************************************
    //  Import model datas
    //*****************************************************

	//----------Model Path----------
	std::string model_path = "sample/navierstokes/model2/";
	
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

    //----------Get incremental Dirichlet condition----------
    int incmax = 10;
    std::vector<double> dufixed = std::vector<double>(ufixed.size());
    for(int i = 0; i < ufixed.size(); i++){
        dufixed[i] = ufixed[i]/(double)incmax;
    }
    
    //----------Initialise with result----------
    std::vector<Vector<double>> u = std::vector<Vector<double>>(nodes.size(), Vector<double>(2));
	std::vector<double> p = std::vector<double>(nodes.size(), 0.0);
    

    //*****************************************************
    //  Update result with NavierStokes equation
    //*****************************************************
    
    //----------Incremental Step----------
    for (int k = 0; k < 1; k++) {
        //----------Culculate Ke Fe and Assembling----------
        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);			//System stiffness matrix
        std::vector<double> Q = std::vector<double>(KDEGREE, 0.0);		//Interna load vector
        for (int i = 0; i < elementsu.size(); i++) {
            Matrix<double> Ke;
            Vector<double> Qe;
            NavierStokes<double, ShapeFunction6Triangle, ShapeFunction3Triangle, Gauss3Triangle>(Ke, Qe, nodes, u, p, elementsu[i], elementsp[i], 1.0);
            Assembling(K, Q, Ke, Qe, elementsu[i], field);
        }

        //----------Culculate traction term and Assembling----------
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);		//External load vector
        /*for(int i = 0; i < edgesu.size(); i++){
            Vector<double> Fe;
            StokesTraction<double, ShapeFunction3Line, ShapeFunction2Line, Gauss2Line>(Fe, nodes, edgesu[i], edgesp[i], 1.0/3.0, 0.0);
            Assembling(F, Fe, edgesu[i], field);
        }*/

        //----------Get residual load vector----------
        std::vector<double> R = std::vector<double>(KDEGREE, 0.0);		//Residual load vector
        for(int i = 0; i < KDEGREE; i++){
            R[i] = F[i] - Q[i];
        }
                    
        //----------Set Dirichlet Boundary Condition----------
        SetDirichlet(K, R, isufixed, dufixed, 1.0e9);

        std::cout << K << std::endl;
        for(auto ri : R){
            std::cout << ri << std::endl;
        }

        //----------Check convergence----------
        double normR = std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
        //double normF = std::inner_product(F.begin(), F.end(), F.begin(), 0.0);
        //std::cout << "k = " << k << "\tNorm = " << normR << std::endl;
        if (normR < 1.0e-6) {
            //std::cout << "k = " << k << "\tNorm = " << normR << std::endl;
            //break;
        }
        
        //----------Solve System Equation----------
        CSR<double> Kmod = CSR<double>(K);
        std::vector<double> result = BiCGSTAB(Kmod, R, 1, 1.0e-20);

        //----------Update displacement u----------
        std::vector<Vector<double> > dup = std::vector<Vector<double> >(nodes.size());
        FieldResultToNodeValue(result, dup, field);
        for(int i = 0; i < nodes.size(); i++){
            u[i] += dup[i].Segment(0, 2);
            if(dup[i].SIZE() == 3){
                p[i] += dup[i](2);
            }
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