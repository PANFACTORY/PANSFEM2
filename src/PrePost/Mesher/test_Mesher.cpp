#include <iostream>
#include <vector>


#include "SquareMesh.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    SquareMesh<double> mesh = SquareMesh<double>(5.0, 4.0, 5, 4);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::vector<int> > edges = mesh.GenerateEdges();
    /*std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(_x.Norm() < 1.0 + 1.0e-5) {
            return true;
        }
        return false;
    });*/
    std::vector<int> elements0 = mesh.GenerateElementIdsSelected([](Vector<double> _x){
        double eps = 1.0e-5;
        if(2.0 - eps < _x(0) && _x(0) < 4.0 + eps && 1.0 - eps < _x(1) && _x(1) < 3.0 + eps) {
            return true;
        }
        return false;
    });

    for(auto i : elements0) {
        std::cout << i << std::endl;    
    }

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();
    return 0;
}