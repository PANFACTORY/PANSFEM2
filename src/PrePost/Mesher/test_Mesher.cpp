#include <iostream>
#include <vector>


#include "SquareAnnulusMesh.h"
#include "../Export/ExportToVTK.h"
#include "../../FEM/Equation/General.h"
#include "../../FEM/Controller/GaussIntegration.h"
#include "../../FEM/Controller/ShapeFunction.h"


using namespace PANSFEM2;


int main(){
    SquareAnnulusMesh2<double> mesh = SquareAnnulusMesh2<double>(10.0, 7.0, 10, 7, 4, 3);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    //std::vector<std::vector<int> > edges = mesh.GenerateEdges();

    double area = 0.0;
    for(auto element : elements) {
        area += Area<double, ShapeFunction4Square, Gauss4Square>(nodes, element);
    }
    std::cout << area << std::endl;

    /*for(auto element : elements) {
        for(auto nodeid : element) {
            std::cout << nodeid << "\t";
        }
        std::cout << std::endl;
    }*/

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();
    return 0;
}