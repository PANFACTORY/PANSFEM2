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

	//----------Initialize design variables----------
	std::vector<double> s = std::vector<double>(elements.size(), 0.5);

	//----------Define parameters----------
	double E0 = 0.001;
	double E1 = 210000.0;				
	double Poisson = 0.3;
	double rho0 = 0.0;
	double rho1 = 1.0;			
	double Density = 0.0000078;	
	double p = 3.0;			

	double dt = 0.001;					
	int step = 101;						

	double beta = 0.3333333333;			
	double ganma = 0.5;					

	double iota = 1.0;
	double lambdamin = 1.0e-20;
	double lambdamax = 1.0e20;
	double lambdaeps = 1.0e-15;
	double movelimit = 0.05;

	double weightlimit = 0.5;
	double objectivebefore = 0.0;
	double objectiveeps = 1.0e-5;

	//----------Optimize loop----------
	for(int k = 0; k < 100; k++){
		std::cout << "k = " << k << "\t";


		//*************************************************
		//  Get weight value and sensitivities
		//*************************************************
        double weight = 0.0;													//Function values of weight
		std::vector<double> dweights = std::vector<double>(s.size());			//Sensitivities of weight
        for (int i = 0; i < elements.size(); i++) {						
			weight += rho0*(1.0 - s[i]) + rho1*s[i] - weightlimit;
			dweights[i] = - rho0 + rho1;
		}


		//*************************************************
		//  Direct time step with Newmark beta method
		//*************************************************
		std::cout << "--->\t"; 
		std::vector<std::vector<Vector<double> > > d = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Displacement of nodes
		std::vector<std::vector<Vector<double> > > v = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Velocity of nodes
		std::vector<std::vector<Vector<double> > > a = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Acceraration of nodes

		for(int t = 1; t < step; t++){
			//----------Assembling----------
			LILCSR<double> A = LILCSR<double>(KDEGREE, KDEGREE);
			std::vector<double> b = std::vector<double>(KDEGREE, 0.0);
			for (int i = 0; i < elements.size(); i++) {
				Vector<double> dne = Vector<double>();
				Vector<double> vne = Vector<double>();
				Vector<double> ane = Vector<double>();
				for(int j = 0; j < elements[i].size(); j++){
					dne = dne.Vstack(d[t - 1][elements[i][j]]);
					vne = vne.Vstack(v[t - 1][elements[i][j]]);
					ane = ane.Vstack(a[t - 1][elements[i][j]]);
				}

				Matrix<double> Ke, Me;
				double E = E0*(1.0 - pow(s[i], p)) + E1*pow(s[i], p);
				double D = Density*(rho0*(1.0 - pow(s[i], p)) + rho1*pow(s[i], p));
				PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
				PlaneMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodes, elements[i], D, 1.0);

				Matrix<double> Ae = Me + pow(dt, 2.0)*beta*Ke;
				Vector<double> be = -Ke*(dne + dt*vne + pow(dt, 2.0)*(0.5 - beta)*ane);

				Assembling(A, b, Ae, be, elements[i], field);
			}

			//----------Set Neumann Boundary Condition----------
			if(t == 1){
				SetNeumann(b, isqfixed, qfixed);
			}
			
			//----------Set Dirichlet Boundary Condition----------
			SetDirichlet(A, b, isufixed, ufixed, 1.0e10);
		
			//----------Solve linear system----------
			CSR<double> Amod = CSR<double>(A);	
			std::vector<double> results = ScalingCG(Amod, b, 100000, 1.0e-10);

			//----------Get d, v, a at step n+1----------
			FieldResultToNodeValue(results, a[t], field);
			for(int i = 0; i < nodes.size(); i++){
				d[t][i] = d[t - 1][i] + dt*v[t - 1][i] + pow(dt, 2.0)*(0.5 - beta)*a[t - 1][i] + pow(dt, 2.0)*beta*a[t][i];
				v[t][i] = v[t - 1][i] + dt*(1.0 - ganma)*a[t - 1][i] + dt*ganma*a[t][i];
			}			
		}


		//*************************************************
		//  Invert time step with Newmark beta method
		//*************************************************
		std::cout << "<---\t"; 
		std::vector<std::vector<Vector<double> > > l = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Displacement of nodes
		std::vector<std::vector<Vector<double> > > m = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Velocity of nodes
		std::vector<std::vector<Vector<double> > > n = std::vector<std::vector<Vector<double> > >(step, std::vector<Vector<double> >(nodes.size(), Vector<double>(2)));	//	Acceraration of nodes

		for(int t = 1; t < step; t++){
			//----------Assembling----------
			LILCSR<double> A = LILCSR<double>(KDEGREE, KDEGREE);
			std::vector<double> b = std::vector<double>(KDEGREE, 0.0);
			for (int i = 0; i < elements.size(); i++) {
				Vector<double> lne = Vector<double>();
				Vector<double> mne = Vector<double>();
				Vector<double> nne = Vector<double>();
				Vector<double> dne = Vector<double>();
 				for(int j = 0; j < elements[i].size(); j++){
					lne = lne.Vstack(l[t - 1][elements[i][j]]);
					mne = mne.Vstack(m[t - 1][elements[i][j]]);
					nne = nne.Vstack(n[t - 1][elements[i][j]]);
					dne = dne.Vstack(d[step - t][elements[i][j]]);
				}

				Matrix<double> Ke, Me;
				double E = E0*(1.0 - pow(s[i], p)) + E1*pow(s[i], p);
				double D = Density*(rho0*(1.0 - pow(s[i], p)) + rho1*pow(s[i], p));
				PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], E, Poisson, 1.0);
				PlaneMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodes, elements[i], D, 1.0);

				Matrix<double> Ae = Me + pow(dt, 2.0)*beta*Ke;
				Vector<double> be = - Ke*dne - Ke*(lne + dt*mne + pow(dt, 2.0)*(0.5 - beta)*nne);

				Assembling(A, b, Ae, be, elements[i], field);
			}

			//----------Set Dirichlet Boundary Condition----------
			SetDirichlet(A, b, isufixed, ufixed, 1.0e10);
		
			//----------Solve linear system----------
			CSR<double> Amod = CSR<double>(A);	
			std::vector<double> results = ScalingCG(Amod, b, 100000, 1.0e-10);

			//----------Get d, v, a at step n+1----------
			FieldResultToNodeValue(results, n[t], field);
			for(int i = 0; i < nodes.size(); i++){
				l[t][i] = l[t - 1][i] + dt*m[t - 1][i] + pow(dt, 2.0)*(0.5 - beta)*n[t - 1][i] + pow(dt, 2.0)*beta*n[t][i];
				m[t][i] = m[t - 1][i] + dt*(1.0 - ganma)*n[t - 1][i] + dt*ganma*n[t][i];
			}			
		}


		//*************************************************
		//	Integration by time
		//*************************************************
		std::cout << "<--->\t"; 
		double objective = 0.0;													//Function value of compliance
		std::vector<double> dobjectives = std::vector<double>(s.size(), 0.0);	//Sensitivities of compliance
		for(int t = 0; t < step; t++){
			for (int i = 0; i < elements.size(); i++) {
				Vector<double> dne = Vector<double>();
				Vector<double> ane = Vector<double>();
				Vector<double> lne = Vector<double>();
				for(int j = 0; j < elements[i].size(); j++){
					dne = dne.Vstack(d[t][elements[i][j]]);
					ane = ane.Vstack(a[t][elements[i][j]]);
					lne = lne.Vstack(l[step - t - 1][elements[i][j]]);
				}

				Matrix<double> Ke, Me;
				double E = E0*(1.0 - pow(s[i], p)) + E1*pow(s[i], p);
				double dE = p*(- E0 + E1)*pow(s[i], p - 1.0);
				double dD = Density*p*(- rho0 + rho1)*pow(s[i], p - 1.0);
				PlaneStrain<double, ShapeFunction8Square, Gauss9Square >(Ke, nodes, elements[i], 1.0, Poisson, 1.0);
				PlaneMass<double, ShapeFunction8Square, Gauss9Square >(Me, nodes, elements[i], 1.0, 1.0);

				if(t == 0 || t == step - 1){
					objective += 0.5*0.5*dne*(E*Ke*dne)*dt;
					dobjectives[i] += 0.5*0.5*dne*(dE*Ke*dne)*dt + 0.5*lne*(dD*Me*ane + dE*Ke*dne)*dt;
				} else {
					objective += 0.5*dne*(E*Ke*dne)*dt;
					dobjectives[i] += 0.5*dne*(dE*Ke*dne)*dt + lne*(dD*Me*ane + dE*Ke*dne)*dt;
				}
			}
		}


		//*************************************************
		//  Post Process
		//*************************************************
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 23), fout);
		AddElementScalers(s, "s", fout, true);
		fout.close();


		//*************************************************
        //  Update design variables with OC method
        //*************************************************

		//----------Check convergence----------
        std::cout << "Objective:\t" << objective << "\tWeight:\t" << weight << "\t";
		if(fabs((objective - objectivebefore) / (objective + objectivebefore)) < objectiveeps) {
			std::cout << std::endl << "----------Convergence----------" << std::endl;
			break;
		}
		
		//----------Get updated design variables with OC method----------
		double lambda0 = lambdamin, lambda1 = lambdamax, lambda;
		std::vector<double> snext = std::vector<double>(elements.size());
		while((lambda1 - lambda0) / (lambda1 + lambda0) > lambdaeps){
			lambda = 0.5 * (lambda1 + lambda0);

			for (int i = 0; i < s.size(); i++) {
				snext[i] = pow(dobjectives[i] / (dweights[i] * lambda), iota) * s[i];
				if(snext[i] < std::max(0.0, (1.0 - movelimit)*s[i])) {
					snext[i] = std::max(0.0, (1.0 - movelimit)*s[i]);
				} else if(snext[i] > std::min(1.0, (1.0 + movelimit)*s[i])) {
					snext[i] = std::min(1.0, (1.0 + movelimit)*s[i]);
				}
			}

			double weightnext = 0.0;
			for (int i = 0; i < elements.size(); i++) {
				weightnext += rho0*(1.0 - snext[i]) + rho1*snext[i] - weightlimit;
			}

			if (weightnext > 0.0) {
				lambda0 = lambda;
			}
			else {
				lambda1 = lambda;
			}
		}

		std::cout << "Lagrange value:\t" << lambda << std::endl;

		//----------Update design variables and objective----------
		for (int i = 0; i < s.size(); i++) {
			s[i] = snext[i];
		}
		objectivebefore = objective;
	}

	
	
	return 0;
}