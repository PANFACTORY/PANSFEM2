#include <iostream>
#include <vector>


#include "Delaunay.h"
#include "../Export/ExportToVTK.h"


using namespace PANSFEM2;


int main(){
    std::vector<Vector<double> > x = { { 0, 0.5 }, { 1, 0.5 }, { 1, 0 }, { 4, 0 }, { 4, 1 }, { 0, 1 } };
    Delaunay<double> mesher = Delaunay<double>(x, { { 5, 4, 3, 2, 1, 0 } }, {}, 0.4);
    std::vector<Vector<double> > nodes = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesher.GenerateFixedlist({ 0 }, [](Vector<double> _x) {
        if(_x(0) < 1.0e-5) {
            return true;
        }
        return false;
    });

    for(auto ui : ufixed) {
        std::cout << ui.first.first << std::endl;
    }

    std::ofstream fout("result.vtk");
    MakeHeadderToVTK(fout);
    AddPointsToVTK(nodes, fout);
    AddElementToVTK(elements, fout);
    AddElementTypes(std::vector<int>(elements.size(), 5), fout);
    fout.close();

    return 0;
}