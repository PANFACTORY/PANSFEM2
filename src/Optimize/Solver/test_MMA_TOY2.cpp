//*********************************************************
//  Minimize    x0^2 + x1^2 + x2^2
//  Subject to  (x0 - 5)^2 + (x1 - 2)^2 + (x2 - 1)^2 <= 9
//              (x0 - 3)^2 + (x1 - 4)^2 + (x2 - 3)^2 <= 9
//              0 <= xj <= 5
//*********************************************************


#include <iostream>
#include <vector>
#include <cmath>


#include "MMA.h"


using namespace PANSFEM2;


int main() {
	//----------Initialize optimization solver----------
	std::vector<double> s = std::vector<double>({4.0, 3.0, 2.0});
	MMA<double> optimizer = MMA<double>(3, 2, 
		1.0, 
		std::vector<double>(2, 0.0), 
		std::vector<double>(2, 1000.0), 
		std::vector<double>(2, 1.0), 
		std::vector<double>(3, 0.0), 
		std::vector<double>(3, 5.0)
	);
	
	//----------Optimize loop----------
	for(int k = 0; k < 10; k++){
		std::cout << std::endl << "k = " << k << "\t";

		//**************************************************
		//	Update design variables
		//**************************************************
		double objective;											                                        //Function value of compliance
		std::vector<double> dobjective = std::vector<double>(3);               								//Sensitivities of compliance
		std::vector<double> constraints = std::vector<double>(2);											//Function values of weight
		std::vector<std::vector<double> > dconstraints = std::vector<std::vector<double> >(2, std::vector<double>(3));	    //Sensitivities of weight

		//----------Get function values and sensitivities----------
		//  Objective
        objective = pow(s[0], 2.0) + pow(s[1], 2.0) + pow(s[2], 2.0); 
        dobjective[0] = 2.0*s[0];
        dobjective[1] = 2.0*s[1];
        dobjective[2] = 2.0*s[2];
        

        //  Constraint1
        constraints[0] = pow(s[0] - 5.0, 2.0) + pow(s[1] - 2.0, 2.0) + pow(s[2] - 1.0, 2.0) - 9.0;
        dconstraints[0][0] = 2.0*(s[0] - 5.0);
        dconstraints[0][1] = 2.0*(s[1] - 2.0);
        dconstraints[0][2] = 2.0*(s[2] - 1.0);


        //  Constraint2
        constraints[1] = pow(s[0] - 3.0, 2.0) + pow(s[1] - 4.0, 2.0) + pow(s[2] - 3.0, 2.0) - 9.0;
        dconstraints[1][0] = 2.0*(s[0] - 3.0);
        dconstraints[1][1] = 2.0*(s[1] - 4.0);
        dconstraints[1][2] = 2.0*(s[2] - 3.0);
        

		std::cout << "f:\t" << objective << "\t";
		std::cout << "g1:\t" << constraints[0] << "\t";
        std::cout << "g2:\t" << constraints[1] << "\t";

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