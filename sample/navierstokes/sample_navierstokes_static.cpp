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

#include "../../src/FEM/Equation/Geometric.h"


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

	//----------Add Field for velocity----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");


    //*****************************************************
    //  Initialize result with Stokes
    //*****************************************************

	//----------Culculate Ke Fe and Assembling----------
    LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);        //  System stiffness matrix
    std::vector<double> F = std::vector<double>(KDEGREE, 0.0);  //  External load vector
    for (int i = 0; i < elementsu.size(); i++) {
        Matrix<double> Ke;
        Stokes<double, ShapeFunction6Triangle, ShapeFunction3Triangle, Gauss3Triangle>(Ke, nodes, elementsu[i], elementsp[i], 1.0);
        Assembling(K, Ke, elementsu[i], field);
    }

    //----------Set Boundary Condition----------
    SetDirichlet(K, F, isufixed, ufixed, 1.0e9);

    //----------Solve System Equation----------
    CSR<double> Kmod = CSR<double>(K);
    std::vector<double> result = CG(Kmod, F, 100000, 1.0e-10);
    std::vector<Vector<double> > up = std::vector<Vector<double> >(nodes.size());     //  Velocity
    FieldResultToNodeValue(result, up, field);
	std::vector<Vector<double>> u = std::vector<Vector<double>>(nodes.size());
	std::vector<double> p = std::vector<double>(nodes.size(), 0.0);
	for(int i = 0; i < nodes.size(); i++){
		u[i] = up[i].Segment(0, 2);
		if(up[i].SIZE() == 3){
			p[i] = up[i](2);
		}
	}


    //*****************************************************
    //  Update result with NavierStokes equation
    //*****************************************************
    
    //----------Incremental Step----------
    for (int k = 0; k < 10; k++) {
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
        
        //----------Get residual load vector----------
        std::vector<double> R = std::vector<double>(KDEGREE, 0.0);		//Residual load vector
        for(int i = 0; i < KDEGREE; i++){
            R[i] = F[i] - Q[i];
        }
                    
        //----------Set Dirichlet Boundary Condition----------
        std::vector<double> tmpufixed = std::vector<double>(isufixed.size(), 0.0);
        SetDirichlet(K, R, isufixed, tmpufixed, 1.0e9);

        //----------Check convergence----------
        double normR = std::inner_product(R.begin(), R.end(), R.begin(), 0.0);
        double normQ = std::inner_product(Q.begin(), Q.end(), Q.begin(), 0.0);
        std::cout << "k = " << k << "\tR Norm = " << normR << "\tQ Norm = " << normQ << std::endl;
        if (normR < 1.0e-3) {
            //std::cout << "k = " << k << "\tNorm = " << normR << std::endl;
            break;
        }
        
        //----------Solve System Equation----------
        CSR<double> Kmod = CSR<double>(K);
        std::vector<double> result = BiCGSTAB(Kmod, R, 10000, 1.0e-10);

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