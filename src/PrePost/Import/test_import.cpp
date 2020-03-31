#include "ImportFromVTK.h"


using namespace PANSFEM2;


int main(){
    ImportModelFromVTK<double> modelimport("../../../sample/planestrain/result.vtk");
    std::vector<Vector<double> > nodes = modelimport.ImportPOINTS();
    std::vector<std::vector<int> > elements = modelimport.ImportCELLS();
    std::vector<double> s = modelimport.ImportCELLSCALARS("s");
    std::vector<Vector<double> > u = modelimport.ImportPOINTVECTORS("u");

    for(auto node : nodes){
        std::cout << node.Transpose();
    }

    std::cout << std::endl;

    for(auto element : elements){
        for(auto node : element){
            std::cout << node << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(auto ui : u){
        std::cout << ui.Transpose();
    }

    return 0;
}