#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/LinearAlgebra/Models/LILCSR.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/FEM/Equation/Fluid.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/IntegrationConstant.h"


using namespace PANSFEM2;


int main() {
	//----------Model Path----------
	std::string model_path = "sample/fluid/";
	
	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
		
	//----------Add Elements for velocity----------
	std::vector<std::vector<int> > elementsu;
	ImportElementsFromCSV(elementsu, model_path + "ElementU.csv");

	//----------Add Elements for pressure----------
	std::vector<std::vector<int> > elementsp;
	ImportElementsFromCSV(elementsp, model_path + "ElementP.csv");

	int elementssize = elementsu.size();
	assert(elementssize == elementsp.size());

	//----------Add Field for velocity----------
	std::vector<int> fieldu;
	int KDEGREEU = 0;
	ImportFieldFromCSV(fieldu, KDEGREEU, nodes.size(), model_path + "FieldU.csv");

    //----------Add Field for pressure----------
	std::vector<int> fieldp;
	int KDEGREEP = 0;
	ImportFieldFromCSV(fieldp, KDEGREEP, nodes.size(), model_path + "FieldP.csv");
	
	//----------Add Dirichlet Condition for velocity----------
	std::vector<int> isufixedu;
	std::vector<double> ufixedu;
	ImportDirichletFromCSV(isufixedu, ufixedu, fieldu, model_path + "DirichletU.csv");

    //----------Add Dirichlet Condition for pressure----------
	std::vector<int> isufixedp;
	std::vector<double> ufixedp;
	ImportDirichletFromCSV(isufixedp, ufixedp, fieldp, model_path + "DirichletP.csv");

	//----------Initialize velocity and pressure----------
	std::vector<Vector<double> > u = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));     //  Velocity
    std::vector<Vector<double> > v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));     //  Auxiliary velocity
    std::vector<double> p = std::vector<double>(nodes.size(), 0.0);                         			//  Pressure
	ImportInitialFromCSV(u, model_path + "DirichletU.csv");

	//----------Define parameters----------
	double dt = 0.01;	//  Time step
	double Re = 10.0;	//	Reynolds number

	//----------Time step loop----------
	for(int t = 0; t < 100; t++){
		std::cout << "t = " << t << std::endl;


        //*************************************************
        //  Get auxiliary velocity
        //*************************************************
		std::cout << "\tGet auxiliary velocity...";

		//----------Culculate Ke Fe and Assembling----------
		LILCSR<double> KV = LILCSR<double>(KDEGREEU, KDEGREEU);         //  System stiffness matrix
		std::vector<double> FV = std::vector<double>(KDEGREEU, 0.0);	//  External load vector
		for (int j = 0; j < elementssize; j++) {
			Vector<double> ue = Vector<double>();
			for(int i = 0; i < elementsu[j].size(); i++){
				ue = ue.Vstack(u[i]);
			}

			Matrix<double> MassTerm;
			FluidMass<double, ShapeFunction8Square, Gauss9Square>(MassTerm, nodes, elementsu[j]);
			
            Matrix<double> ConvectiveTerm;
            Convective<double, ShapeFunction8Square, Gauss9Square>(ConvectiveTerm, nodes, u, elementsu[j]);

            Matrix<double> DiffusionTerm;
            Diffusion<double, ShapeFunction8Square, Gauss9Square>(DiffusionTerm, nodes, elementsu[j]);

            Vector<double> Fe = (MassTerm - dt*(ConvectiveTerm + DiffusionTerm/Re))*ue;

			Assembling(KV, FV, MassTerm, Fe, elementsu[j], fieldu);
		}

		//----------Set Boundary Condition----------
		SetDirichlet(KV, FV, isufixedu, ufixedu, 1.0e10);

		//----------Solve System Equation----------
		CSR<double> KVmod = CSR<double>(KV);
		std::vector<double> resultv = ScalingCG(KVmod, FV, 100000, 1.0e-10);
		FieldResultToNodeValue(resultv, v, fieldu);

		std::cout << "done" << std::endl;

        //*************************************************
        //  Solve Poisson equation for pressure
        //*************************************************
		std::cout << "\tGet pressure distribution...";

        //----------Culculate Ke Fe and Assembling----------
		LILCSR<double> KP = LILCSR<double>(KDEGREEP, KDEGREEP);         //  System stiffness matrix
		std::vector<double> FP = std::vector<double>(KDEGREEP, 0.0);	//  External load vector
		for (int j = 0; j < elementssize; j++) {
			Matrix<double> Ke;
            Vector<double> Fe;
			Poisson<double, ShapeFunction8Square, ShapeFunction4Square, Gauss4Square>(Ke, Fe, nodes, v, elementsu[j], elementsp[j]);
			Fe /= -dt;
			Assembling(KP, FP, Ke, Fe, elementsp[j], fieldp);
		}

		//----------Set Boundary Condition----------
		SetDirichlet(KP, FP, isufixedp, ufixedp, 1.0e10);

		//----------Solve System Equation----------
		CSR<double> KPmod = CSR<double>(KP);
		std::vector<double> resultp = ScalingCG(KPmod, FP, 100000, 1.0e-10);
		FieldResultToNodeValue(resultp, p, fieldp);

		std::cout << "done" << std::endl;

        //*************************************************
        //  Get velocities of next step
        //*************************************************
        std::cout << "\tGet velocity...";

		//----------Culculate Ke Fe and Assembling----------
		LILCSR<double> KU = LILCSR<double>(KDEGREEU, KDEGREEU);         //  System stiffness matrix
		std::vector<double> FU = std::vector<double>(KDEGREEU, 0.0);	//  External load vector
		for (int j = 0; j < elementssize; j++) {
			Vector<double> pe = Vector<double>(elementsp[j].size());
			for(int i = 0; i < elementsp[j].size(); i++){
				pe(i) = p[elementsp[j][i]];
			}

            Vector<double> ve = Vector<double>();
			for(int i = 0; i < elementsu[j].size(); i++){
				ve = ve.Vstack(v[elementsu[j][i]]);
			}
			
			Matrix<double> Me;
			FluidMass<double, ShapeFunction8Square, Gauss9Square>(Me, nodes, elementsu[j]);
			
            Matrix<double> Ke;
            UpdateVelocity<double, ShapeFunction8Square, ShapeFunction4Square, Gauss9Square>(Ke, nodes, elementsu[j], elementsp[j]);

            Vector<double> Fe = Me*ve - dt*Ke*pe;

			Assembling(KU, FU, Me, Fe, elementsu[j], fieldu);
		}

		//----------Set Boundary Condition----------
		SetDirichlet(KU, FU, isufixedu, ufixedu, 1.0e10);

		//----------Solve System Equation----------
		CSR<double> KUmod = CSR<double>(KU);
		std::vector<double> resultu = ScalingCG(KUmod, FU, 100000, 1.0e-10);
		FieldResultToNodeValue(resultu, u, fieldu);

		std::cout << "done" << std::endl;

        //*************************************************		
		//  Save values
        //*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elementsp, fout);
		AddElementTypes(std::vector<int>(elementsp.size(), 9), fout);
        AddPointVectors(u, "u", fout, true);
		AddPointVectors(v, "v", fout, false);
		AddPointScalers(p, "p", fout, false);
		AddPointScalers(FP, "fp", fout, false);
		fout.close();
	}

	return 0;
}