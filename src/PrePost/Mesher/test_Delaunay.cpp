#include <iostream>
#include <vector>


#include "Delaunay.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    std::vector<Vector<double> > x = { { 0, 0.5 }, { 1, 0.5 }, { 1, 0 }, { 4, 0 }, { 4, 1 }, { 0, 1 } };
    Delaunay<double> mesher = Delaunay<double>(x, { { 5, 4, 3, 2, 1, 0 } }, {}, 1.0);
    std::vector<Vector<double> > nodes = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();

    return 0;
}