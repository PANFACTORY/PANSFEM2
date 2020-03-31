#include <iostream>
#include <vector>


#include "StructuredGrid.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    StructuredGrid<double> mesh = StructuredGrid<double>(4, 3, 4, 3);
    std::vector<Vector<double> > nodes = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();
    return 0;
}