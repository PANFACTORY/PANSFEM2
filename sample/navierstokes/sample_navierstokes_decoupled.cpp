#include <iostream>
#include <vector>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/NavierStokes.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"
#include "../../src/FEM/Equation/General.h"


using namespace PANSFEM2;


int main() {
    //----------パラメータ設定----------
    double rho = 1.0;
    double mu = 1.0/100.0;
    int tmax = 800;
    double dt = 0.005;
    

    //----------モデリング----------
    SquareMesh<double> mesher = SquareMesh<double>(1, 1, 30, 30);
    std::vector<Vector<double> > x = mesher.GenerateNodes();
    std::vector<std::vector<int> > elements = mesher.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed0 = mesher.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x) {
        if(fabs(_x(0) - 0.0) < 1.0e-5 || fabs(_x(1) - 0.0) < 1.0e-5 || fabs(_x(0) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > ufixed1 = mesher.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x) {
        if(fabs(_x(1) - 1.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& ufixed1i : ufixed1) {
        if(ufixed1i.first.second == 0) {
            ufixed1i.second = 1.0;
        }
    }
    std::vector<std::pair<std::pair<int, int>, double> > pfixed = mesher.GenerateFixedlist({ 0 }, [](Vector<double> _x) {
        if(fabs(_x(0) - 0.5) < 1.0e-5 && fabs(_x(1)) < 1.0e-5) {
            return true;
        }
        return false;
    });


    //----------変数定義と境界条件適用----------
    std::vector<Vector<double> > v = std::vector<Vector<double> >(x.size(), Vector<double>(2));
    std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
	std::vector<std::vector<int> > nodetoglobalvu = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
    SetDirichlet(v, nodetoglobalvu, ufixed0);
    SetDirichlet(u, nodetoglobalvu, ufixed0);
    SetDirichlet(v, nodetoglobalvu, ufixed1);
    SetDirichlet(u, nodetoglobalvu, ufixed1);
	int KDEGREEVU = Renumbering(nodetoglobalvu);

    std::vector<Vector<double> > p = std::vector<Vector<double> >(x.size(), Vector<double>(1));
    std::vector<std::vector<int> > nodetoglobalp = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
    SetDirichlet(p, nodetoglobalp, pfixed);
	int KDEGREEP = Renumbering(nodetoglobalp);


    //----------Time step loop----------
    for(int t = 0; t < tmax; t++) {
        std::cout << "t=" << t << std::endl;

        //----------Get auxiliary velocity----------
        LILCSR<double> Kv = LILCSR<double>(KDEGREEVU, KDEGREEVU);		//  System stiffness matrix
        std::vector<double> Fv = std::vector<double>(KDEGREEVU, 0.0);	//  System load vector

        for (int i = 0; i < elements.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Vector<double> Fe;
            NavierStokesAuxiliaryVelocity<double, ShapeFunction4Square, Gauss4Square>(Ke, Fe, nodetoelement, elements[i], { 0, 1 }, x, u, rho, mu, dt);
            Assembling(Kv, Fv, v, Ke, nodetoglobalvu, nodetoelement, elements[i]);
            Assembling(Fv, Fe, nodetoglobalvu, nodetoelement, elements[i]);
        }

        CSR<double> Kvmod = CSR<double>(Kv);
        std::vector<double> resultv = ScalingCG(Kvmod, Fv, 100000, 1.0e-10);
        Disassembling(v, resultv, nodetoglobalvu);


        //----------Solve Pressure-Poisson equation----------
        LILCSR<double> Kp = LILCSR<double>(KDEGREEP, KDEGREEP);			//  System stiffness matrix
        std::vector<double> Fp = std::vector<double>(KDEGREEP, 0.0);	//  System load vector

        for (int i = 0; i < elements.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Vector<double> Fe;
            NavierStokesPressurePoisson<double, ShapeFunction4Square, Gauss4Square>(Ke, Fe, nodetoelement, elements[i], { 0 }, x, v, rho, dt);
            Assembling(Kp, Fp, p, Ke, nodetoglobalp, nodetoelement, elements[i]);
            Assembling(Fp, Fe, nodetoglobalp, nodetoelement, elements[i]);
        }

        CSR<double> Kpmod = CSR<double>(Kp);
        std::vector<double> resultp = ScalingCG(Kpmod, Fp, 100000, 1.0e-10);
        Disassembling(p, resultp, nodetoglobalp);


        //----------Update next step velocity----------
        LILCSR<double> Ku = LILCSR<double>(KDEGREEVU, KDEGREEVU);		//  System stiffness matrix
        std::vector<double> Fu = std::vector<double>(KDEGREEVU, 0.0);	//  System load vector

        for (int i = 0; i < elements.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            Vector<double> Fe;
            NavierStokesNextstepVelocity<double, ShapeFunction4Square, Gauss4Square>(Ke, Fe, nodetoelement, elements[i], { 0, 1 }, x, v, p, rho, dt);
            Assembling(Ku, Fu, u, Ke, nodetoglobalvu, nodetoelement, elements[i]);
            Assembling(Fu, Fe, nodetoglobalvu, nodetoelement, elements[i]);
        }

        CSR<double> Kumod = CSR<double>(Ku);
        std::vector<double> resultu = ScalingCG(Kumod, Fu, 100000, 1.0e-10);
        Disassembling(u, resultu, nodetoglobalvu);    
        

        //----------Export result----------
        std::ofstream fout("sample/navierstokes/result" + std::to_string(t) + ".vtk");
        MakeHeadderToVTK(fout);
        AddPointsToVTK(x, fout);
        AddElementToVTK(elements, fout);
        AddElementTypes(std::vector<int>(elements.size(), 5), fout);
        AddPointVectors(u, "u", fout, true);
        AddPointScalers(p, "p", fout, false);
        AddPointVectors(v, "v", fout, false);
        fout.close();
    }
    
	return 0;
}