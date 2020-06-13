#include <iostream>
#include <vector>


#include "GrandStructure.h"
#include "../Export/ExportToVTK.h"
#include "../../FEM/Equation/General.h"
#include "../../FEM/Controller/GaussIntegration.h"
#include "../../FEM/Controller/ShapeFunction.h"


using namespace PANSFEM2;


int main(){
    GrandStructure2D<double> mesh = GrandStructure2D<double>(4.0, 3.0, 4, 3, 1.5);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    //std::vector<std::vector<int> > edges = mesh.GenerateEdges2();
    /*std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(_x.Norm() < 1.0 + 1.0e-5) {
            return true;
        }
        return false;
    });

    /*double area = 0.0;
    for(auto element : elements) {
        area += Area<double, ShapeFunction4Square, Gauss4Square>(nodes, element);
    }
    std::cout << area << std::endl;
    

    for(auto edge : edges) {
        for(auto nodeid : edge) {
            std::cout << nodeid << "\t";
        }
        std::cout << std::endl;    
    }*/

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 3), fout);
    fout.close();
    return 0;
}