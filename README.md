# PANSFEM2
is a FEM library for C++. This contains  
*  Linear algebra models (Vector, Dense Matrix, Sparse Matrix)
*  Linear algebra solvers (CG, BiCGSTAB, BiCGSTAB2, Lanczos, ...)
*  FEM controller (Assembling system, Shape function, Integration method, ...)
*  FEM equations (Linear/Total Lagrange/Updated Lagrange solid, Heat transfer, Fluid, ...)
*  Optimization solver (OC method, MMA)

and you can easily create the simulation environment you need.

## Build and Run
You can build and run a sample code as below

---
$ git clone https://github.com/PANFACTORY/PANSFEM2.git  
$ cd PANSFEM2  
$ g++ -O3 -fopenmp sample/solid/sample_updatedlagrange.cpp  
$ a.exe  

---
and you will get vtk format (*.vtk) result file.  

## Selling point
This library has 3 advantages:

*  Template library  
→ You only pass through the path and include this.
*  Various combination  
→ You can select Equations, Shape functions and Integration methods. 
*  Optimization  
→ This program is originaly developed for Topology Optimization. 

## Samples

![Navier-Stokes](https://github.com/PANFACTORY/PANSFEM2/blob/images/img/NS_Decouple.JPG)

## Author
Tanabe Yuta([@Bena43754971](https://twitter.com/Bena43754971))


## LICENSE
This library is licensed under the [MIT](https://github.com/PANFACTORY/PANSFEM2/blob/master/LICENSE) license.  
(C)2020 TanabeYuta
