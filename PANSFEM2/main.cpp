#include <iostream>
#include <vector>
#include "LinearAlgebra/Models/Point.h"
#include "FEM/Models/Mesh.h"
#include "FEM/Models/Field.h"


using namespace PANSFEM2;


int main() {
	Mesh<double> m = Mesh<double>(2);
	m.setpnode(new Point<double>({ 0.0, 0.0 }));
	m.setpnode(new Point<double>({ 1.0, 0.0 }));
	m.setpnode(new Point<double>({ 1.0, 1.0 }));
	m.setpnode(new Point<double>({ 0.0, 1.0 }));
	m.setelement({ 0, 1, 2, 3 });

	Field<double> field = Field<double>(&m, 1, 1);

	return 0;
}