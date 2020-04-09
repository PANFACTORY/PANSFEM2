#include <iostream>
#include <vector>


#include "Annulus.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    Annulus<double> mesh = Annulus<double>(100.0, 200.0, 100, 720);
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