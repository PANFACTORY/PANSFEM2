#include <iostream>
#include <vector>
#include <cmath>


#include "MMA.h"


using namespace PANSFEM2;


int main() {
    //----------Set constant----------
    double C1 = 1.0;
    double C2 = 0.124;

	//----------Initialize optimization solver----------
	std::vector<double> s = { 1.5, 0.5 };
	MMA<double> optimizer = MMA<double>(2, 2, 1.0, 
		{0.0, 0.0}, 
		{1000.0, 1000.0}, 
		{1.0, 1.0}, 
		{ 0.2, 0.1 }, 
		{ 4.0, 1.6 });
	optimizer.SetParameters(1.0e-5, 0.1, 0.5, 0.5, 0.7, 1.2, 1.0e-5);
	
	//----------Optimize loop----------
	for(int k = 0; k < 100; k++){
		std::cout << std::endl << "k = " << k << "\t";

		//**************************************************
		//	Update design variables
		//**************************************************
		double objective;											                                        	//Function value of compliance
		std::vector<double> dobjective = std::vector<double>(2);               									//Sensitivities of compliance
		std::vector<double> constraints = std::vector<double>(2);												//Function values of weight
		std::vector<std::vector<double> > dconstraints = std::vector<std::vector<double> >(2, std::vector<double>(2));	//Sensitivities of weight

		//----------Get function values and sensitivities----------
		//  Objective
        objective = C1*s[0]*sqrt(1.0 + pow(s[1], 2.0)); 
        dobjective[0] = C1*sqrt(1.0 + pow(s[1], 2.0));
        dobjective[1] = C1*s[0]*s[1]/sqrt(1.0 + pow(s[1], 2.0));

        //  Constraint
        constraints[0] = C2*sqrt(1.0 + pow(s[1], 2.0))*(8.0*pow(s[0], -1.0) + 1.0*pow(s[0], -1.0)*pow(s[1], -1.0)) - 1.0;
        dconstraints[0][0] = C2*sqrt(1.0 + pow(s[1], 2.0))*(-8.0*pow(s[0], -2.0) - 1.0*pow(s[0], -2.0)*pow(s[1], -1.0));
        dconstraints[0][1] = C2*s[1]/sqrt(1.0 + pow(s[1], 2.0))*(8.0*pow(s[0], -1.0) + 1.0*pow(s[0], -1.0)*pow(s[1], -1.0)) + C2*sqrt(1.0 + pow(s[1], 2.0))*(-1.0*pow(s[0], -1.0)*pow(s[1], -2.0));
        constraints[1] = C2*sqrt(1.0 + pow(s[1], 2.0))*(8.0*pow(s[0], -1.0) - 1.0*pow(s[0], -1.0)*pow(s[1], -1.0)) - 1.0;
        dconstraints[1][0] = C2*sqrt(1.0 + pow(s[1], 2.0))*(-8.0*pow(s[0], -2.0) + 1.0*pow(s[0], -2.0)*pow(s[1], -1.0));
        dconstraints[1][1] = C2*s[1]/sqrt(1.0 + pow(s[1], 2.0))*(8.0*pow(s[0], -1.0) - 1.0*pow(s[0], -1.0)*pow(s[1], -1.0)) + C2*sqrt(1.0 + pow(s[1], 2.0))*(1.0*pow(s[0], -1.0)*pow(s[1], -2.0));
        
		std::cout << "Objective:\t" << objective << "\t";
		std::cout << "Constraints0:\t" << constraints[0] << "\t";
        std::cout << "Constraints1:\t" << constraints[1] << "\t";

        for(auto si : s){
            std::cout << "\t" << si; 
        }	

		//----------Check convergence----------
		if(optimizer.IsConvergence(objective)){
			std::cout << std::endl << "--------------------Optimized--------------------" << std::endl;
			break;
		}

		//----------Update s----------
		optimizer.UpdateVariables(s, objective, dobjective, constraints, dconstraints);
	}

	return 0;
}