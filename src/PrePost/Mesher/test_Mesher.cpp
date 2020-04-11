#include <iostream>
#include <vector>


#include "SquareCircleAnnulus.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    SquareCircleAnnulus<double> mesh = SquareCircleAnnulus<double>(1.0, 1.0, 0.3, 8, 8, 4);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<int> field = mesh.GenerateFields(2);
    /*std::vector<int> isufixed = mesh.GenerateFixedlist(2, { 0, 1 }, [](Vector<double> _x){
        if(_x.Norm() < 1.0 + 1.0e-5) {
            return true;
        }
        return false;
    });*/

    /*for(auto i : field) {
        std::cout << i << std::endl;
    }*/

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();
    return 0;
}