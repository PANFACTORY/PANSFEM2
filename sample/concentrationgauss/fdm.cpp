#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/LinearAlgebra/Models/Matrix.h"
#include "../../src/PrePost/Import/ImportFromCSV.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/Optimize/Filter/Heaviside.h"


using namespace PANSFEM2;


int main() {
	//*****************************************************
    //  Import model
    //*****************************************************
	std::string model_path = "sample/concentrationgauss/";
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");
	

    //*****************************************************
	//  Setting filter
    //*****************************************************
    std::vector<Vector<double> > cg = std::vector<Vector<double> >(elements.size());
    for(int i = 0; i < elements.size(); i++){
        Vector<double> centeri = Vector<double>(2);
        for(auto nodei : elements[i]){
            centeri += nodes[nodei];
        }
        cg[i] = centeri/(double)elements[i].size();
    }
    std::vector<std::vector<int> > neighbors = std::vector<std::vector<int> >(elements.size(), std::vector<int>());
    std::vector<std::vector<double> > w = std::vector<std::vector<double> >(elements.size(), std::vector<double>());
    double R = 6.0;
    for(int i = 0; i < elements.size(); i++){
        for(int j = 0; j < elements.size(); j++){
            if((cg[i] - cg[j]).Norm() <= R){
                neighbors[i].push_back(j);
                w[i].push_back((R - (cg[i] - cg[j]).Norm())/R);
            }
        }
    }
    HeavisideFilter<double> filter = HeavisideFilter<double>(elements.size(), neighbors, w);
    filter.UpdateBeta(1.0);


    //*****************************************************
    //  Define parameters
    //*****************************************************
    double voidlimit = 0.2;
    double perimeterlimit = 0.04;
    double perimetereps = 1.0e-10;
    double ds0 = 1.0e-7;
    double ds1 = 1.0e-7;


    //*****************************************************
	//  Initialize design variables
    //*****************************************************
	std::vector<double> s0 = std::vector<double>(elements.size(), 0.5);
    std::vector<double> s1 = std::vector<double>(elements.size(), 0.5);
    std::vector<double> r0 = filter.GetFilteredVariables(s0);
    std::vector<double> r1 = filter.GetFilteredVariables(s1);

    double h = 0.0;
    double l = 0.0;
    for (int i = 0; i < elements.size(); i++) {						
        h += (1.0 - r0[i])/(voidlimit*elements.size());
        for(auto j : neighbors[i]){
            l += (1.0 - r0[j]*r1[j])*sqrt(pow(r0[i] - r0[j], 2.0) + perimetereps)/(double)neighbors[i].size()/(perimeterlimit*elements.size());
        }
    }
    h -= 1.0;
    l -= 1.0;

    std::vector<double> dhdr0 = std::vector<double>(s0.size(), -h);    //Sensitivities of weight
    std::vector<double> dhdr1 = std::vector<double>(s1.size(), -h);    //Sensitivities of weight
    std::vector<double> dldr0 = std::vector<double>(s0.size(), -l);    //Sensitivities of perimeter
    std::vector<double> dldr1 = std::vector<double>(s1.size(), -l);    //Sensitivities of perimeter

   
    //*****************************************************
	//  Get each sensitivity loop
    //*****************************************************
	for(int k = 0; k < elements.size(); k++){
        //----------Get filterd design variables----------
        double tmps0k = s0[k];
        s0[k] += ds0;
        r0 = filter.GetFilteredVariables(s0);
        r1 = filter.GetFilteredVariables(s1);


        //----------Get sensitivities----------
        for (int i = 0; i < elements.size(); i++) {						
            dhdr0[k] += (1.0 - r0[i])/(voidlimit*elements.size());
            for(auto j : neighbors[i]){
                dldr0[k] += (1.0 - r0[j]*r1[j])*sqrt(pow(r0[i] - r0[j], 2.0) + perimetereps)/(double)neighbors[i].size()/(perimeterlimit*elements.size());
            }
		}
        dhdr0[k] -= 1.0;
        dldr0[k] -= 1.0;


        //----------Undo design variables----------
        s0[k] = tmps0k;
        dhdr0[k] /= ds0;
        dldr0[k] /= ds0;
        

        //----------Get filterd design variables----------
        double tmps1k = s1[k];
        s1[k] += ds1;
        r0 = filter.GetFilteredVariables(s0);
        r1 = filter.GetFilteredVariables(s1);


        //----------Get sensitivities----------
        for (int i = 0; i < elements.size(); i++) {						
            dhdr1[k] += (1.0 - r0[i])/(voidlimit*elements.size());
            for(auto j : neighbors[i]){
                dldr1[k] += (1.0 - r0[j]*r1[j])*sqrt(pow(r0[i] - r0[j], 2.0) + perimetereps)/(double)neighbors[i].size()/(perimeterlimit*elements.size());
            }
		}
        dhdr1[k] -= 1.0;
        dldr1[k] -= 1.0;


        //----------Undo design variables----------
        s1[k] = tmps1k;
        dhdr1[k] /= ds1;
        dldr1[k] /= ds1;
	}


    //*************************************************
    //  Filtering sensitivities
    //*************************************************
    std::vector<double> dhds0 = filter.GetFilteredSensitivitis(s0, dhdr0);
    std::vector<double> dhds1 = filter.GetFilteredSensitivitis(s1, dhdr1);
    std::vector<double> dlds0 = filter.GetFilteredSensitivitis(s0, dldr0);
    std::vector<double> dlds1 = filter.GetFilteredSensitivitis(s1, dldr1);


    //*************************************************
    //  Post Process
    //*************************************************
    std::ofstream fout(model_path + "fdm.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 9), fout);
    AddElementScalers(r0, "r0", fout, true);
    AddElementScalers(r1, "r1", fout, false);
    AddElementScalers(dhds0, "dhds0", fout, false);
    AddElementScalers(dhds1, "dhds1", fout, false);
    AddElementScalers(dlds0, "dlds0", fout, false);
    AddElementScalers(dlds1, "dlds1", fout, false);
    fout.close();
	
	return 0;
}