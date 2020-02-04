#include <iostream>
#include <vector>
#include <cmath>


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
	
	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

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

	//----------Add Neumann Condition for velocity----------
	std::vector<int> isqfixedu;
	std::vector<double> qfixedu;
	ImportNeumannFromCSV(isqfixedu, qfixedu, fieldu, model_path + "NeumannU.csv");

    //----------Add Dirichlet Condition for pressure----------
	std::vector<int> isufixedp;
	std::vector<double> ufixedp;
	ImportDirichletFromCSV(isufixedp, ufixedp, fieldp, model_path + "DirichletP.csv");

	//----------Add Neumann Condition for pressure----------
	std::vector<int> isqfixedp;
	std::vector<double> qfixedp;
	ImportNeumannFromCSV(isqfixedp, qfixedp, fieldp, model_path + "NeumannP.csv");

	//----------Initialize velocity and pressure----------
	std::vector<Vector<double> > u = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));     //  Velocity
    std::vector<Vector<double> > v = std::vector<Vector<double> >(nodes.size(), Vector<double>(2));     //  Auxiliary velocity
    std::vector<double> p = std::vector<double>(nodes.size(), 0.0);                         			//  Pressure

	ImportInitialFromCSV(u, model_path + "DirichletU.csv");

	//----------Define parameters----------
	double dt = 0.01;		//  Time step
	double Re = 1.0;		//	Reynolds number

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
		for (auto element : elements) {
			Vector<double> ue = Vector<double>();
			for(auto i : element){
				ue = ue.Vstack(u[i]);
			}

			Matrix<double> MassTerm;
			FluidMass<double, ShapeFunction4Square, Gauss4Square>(MassTerm, nodes, element);
			
            Matrix<double> ConvectiveTerm;
            Convective<double, ShapeFunction4Square, Gauss4Square>(ConvectiveTerm, nodes, u, element);

            Matrix<double> DiffusionTerm;
            Diffusion<double, ShapeFunction4Square, Gauss4Square>(DiffusionTerm, nodes, element);

            Vector<double> Fe = (MassTerm - dt*(ConvectiveTerm + DiffusionTerm/Re))*ue;

			Assembling(KV, FV, MassTerm, Fe, element, fieldu);
		}

		//----------Set Boundary Condition----------
		SetNeumann(FV, isqfixedu, qfixedu);
		SetDirichlet(KV, FV, isufixedu, ufixedu, 1.0e5);

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
		for (auto element : elements) {
			Matrix<double> Ke;
            Vector<double> Fe;
			Poisson<double, ShapeFunction4Square, Gauss4Square>(Ke, Fe, nodes, v, element);
			Fe /= -dt;
			Assembling(KP, FP, Ke, Fe, element, fieldp);
		}

		//----------Set Boundary Condition----------
		SetNeumann(FP, isqfixedp, qfixedp);
		SetDirichlet(KP, FP, isufixedp, ufixedp, 1.0e5);

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
		for (auto element : elements) {
			Vector<double> pe = Vector<double>(element.size());
			for(int i = 0; i < element.size(); i++){
				pe(i) = p[element[i]];
			}

            Vector<double> ve = Vector<double>();
			for(auto i : element){
				ve = ve.Vstack(v[i]);
			}
			
			Matrix<double> Me;
			FluidMass<double, ShapeFunction4Square, Gauss4Square>(Me, nodes, element);
			
            Matrix<double> Ke;
            UpdateVelocity<double, ShapeFunction4Square, Gauss4Square>(Ke, nodes, element);

            Vector<double> Fe = Me*ve - dt*Ke*pe;

			Assembling(KU, FU, Me, Fe, element, fieldu);
		}

		//----------Set Boundary Condition----------
		SetNeumann(FU, isqfixedu, qfixedu);
		SetDirichlet(KU, FU, isufixedu, ufixedu, 1.0e5);

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
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
        AddPointVectors(u, "u", fout, true);
		AddPointVectors(v, "v", fout, false);
		AddPointScalers(p, "p", fout, false);
		fout.close();
	}

	return 0;
}