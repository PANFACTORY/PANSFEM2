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
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/optimize_dynamic/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	
	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Add Neumann Condition----------
	std::vector<int> isqfixed;
	std::vector<double> qfixed;
	ImportNeumannFromCSV(isqfixed, qfixed, field, model_path + "Neumann.csv");

	//----------Add Initial Condition----------
	std::vector<Vector<double> > dn;	//	Displacement of nodes at step n
	std::vector<Vector<double> > vn;	//	Velocity of nodes at step n
	std::vector<Vector<double> > an;	//	Acceraration of nodes at step n
    	
	//----------Define parameters----------
	double E = 210000.0;				//	Young moduls
	double Poisson = 0.3;				//	Poisson ratio
	double rho = 0.0000078;				//	Density

	double dt = 0.01;					//	Time step

	double beta = 0.3333333333;			//	Parameter beta for Newmark beta method
	double ganma = 0.5;					//	Parameter ganma for Newmark beta method
	
	//----------Time step loop----------
	for(int t = 0; t < 100; t++){
		std::cout << "t = " << t << "\t";


        //*************************************************
        //  Direct analysis with Newmark beta method
        //*************************************************
    
        //----------Assembling----------
		LILCSR<double> A = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> b = std::vector<double>(KDEGREE, 0.0);
		for (int i = 0; i < elements.size(); i++) {
			Vector<double> dne = Vector<double>();
			Vector<double> vne = Vector<double>();
            Vector<double> ane = Vector<double>();
            for(int j = 0; j < elements[i].size(); j++){
                dne = dne.Vstack(dn[elements[i][j]]);
                vne = vne.Vstack(vn[elements[i][j]]);
                ane = ane.Vstack(an[elements[i][j]]);
            }

			Matrix<double> Ke, Me;
			PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
			PlaneMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodes, elements[i], rho, 1.0);

			Matrix<double> Ae = Me + pow(dt, 2.0)*beta*Ke;
			Vector<double> be = -Ke*(dne + dt*vne + pow(dt, 2.0)*(0.5 - beta)*ane);

			Assembling(A, b, Ae, be, elements[i], field);
		}

        //----------Set Neumann Boundary Condition----------
	    SetNeumann(b, isqfixed, qfixed);

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(A, b, isufixed, ufixed, 1.0e10);
        
        //----------Solve linear system----------
        CSR<double> Amod = CSR<double>(A);	
        std::vector<double> results = ScalingCG(Amod, b, 100000, 1.0e-10);

		//----------Get d, v, a at step n+1----------
		std::vector<Vector<double> > dnp1;	//	Displacement of nodes at step n+1
		std::vector<Vector<double> > vnp1;	//	Velocity of nodes at step n+1
		std::vector<Vector<double> > anp1;	//	Acceraration of nodes at step n+1

		FieldResultToNodeValue(results, anp1, field);
		for(int i = 0; i < nodes.size(); i++){
			dnp1[i] = dn[i] + dt*vn[i] + pow(dt, 2.0)*(0.5 - beta)*an[i] + pow(dt, 2.0)*beta*anp1[i];
			vnp1[i] = vn[i] + dt*(1.0 - ganma)*an[i] + dt*ganma*anp1[i];
		}


        //*************************************************
        //  Post Process
        //*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 23), fout);
		AddPointVectors(dn, "d", fout, true);
		fout.close();


		//*************************************************
		//	Update d, v, a
		//*************************************************
		dn = dnp1;
		vn = vnp1;
		an = anp1;
	}
	
	return 0;
}