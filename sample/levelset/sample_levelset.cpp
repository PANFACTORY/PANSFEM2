#include <iostream>
#include <vector>
#include <numeric>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/PlaneStrain.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/Optimize/Method/LevelSet.h"


using namespace PANSFEM2;


int main() {
    //----------ステップ0：設計変数の設定----------
    double Vmax = 0.4;
    double tau = 0.07;
    double E0 = 1.0;
    double Emin = 1.0e-4;
    double nu = 0.3;
    double nvol = 100;      //  iteration number for volume constraint
    double dt = 0.1;
    double d = -0.02;
    double p = 4.0;

    double A1 = -1.5*(1.0 - nu)*(1.0 - 14.0*nu + 15.0*pow(nu, 2.0))*E0/((1.0 + nu)*(7.0 - 5.0*nu)*pow(1.0 - 2.0*nu, 2.0));
    double A2 = 7.5*(1.0 - nu)*E0/((1.0 + nu)*(7.0 - 5.0*nu));
    double c = A1/(A1 + 2.0*A2);


    //----------ステップ1：固定設計領域と境界条件の設定----------
    SquareMesh<double> mesh = SquareMesh<double>(80.0, 60.0, 160, 120);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 80.0) < 1.0e-5 && abs(_x(1) - 30.0) < 1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -1.0;
    }


    std::vector<Vector<double> > phi = std::vector<Vector<double> >(x.size(), { 1.0 });     //  φ
    std::vector<double> str = std::vector<double>(elements.size(), 1.0);                    //  χ(φ)
    double volInit = std::accumulate(str.begin(), str.end(), 0.0)/(double)elements.size();
    

    //----------最適化ループ----------
    for(int t = 0; t < 200; t++) {
        //----------ステップ2：固定設計領域の有限要素離散化と数値解析----------
        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        
        SetDirichlet(u, nodetoglobal, ufixed);
        int KDEGREE = Renumbering(nodetoglobal);

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, Emin + str[i]*(E0 - Emin), nu, 1.0);
            Assembling(K, F, u, Ke, nodetoglobal, nodetoelement, elements[i]);
		}

        Assembling(F, qfixed, nodetoglobal);

        CSR<double> Kmod = CSR<double>(K);	
        std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(u, result, nodetoglobal);


        //----------ステップ3：目的汎関数と制約汎関数の計算----------
        std::vector<double> SED = std::vector<double>(elements.size());
        for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, Emin + str[i]*(E0 - Emin), nu, 1.0);
            Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            SED[i] = ue*(Ke*ue);
        }

        double objective = std::accumulate(SED.begin(), SED.end(), 0.0);
        double vol = std::accumulate(str.begin(), str.end(), 0.0)/(double)elements.size();

        std::cout << "t = " << t << "\tCompliance = " << objective << "\tVolume = " << vol << std::endl;


        //----------ステップ4：収束条件の判定----------

        //----------ステップ5：目的汎関数と制約汎関数の設計感度の計算----------
        std::vector<double> TD = std::vector<double>(elements.size());
        for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStrainStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, (A1 + 2.0*A2)*pow(1.0 - c, 2.0), c, 1.0);
            Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            SED[i] = (1.0e-4 + str[i]*(1.0 - 1.0e-4))*ue*(Ke*ue);
        }


        //----------ステップ6：レベルセット関数の更新，ステップ2に戻る----------
        double ex = Vmax + (volInit - Vmax)*std::max(0.0, 1.0 - t/(double)nvol);
        double lambda = std::accumulate(TD.begin(), TD.end(), 0.0)/(double)x.size()*exp(p*((vol - ex)/ex + d));
        double C = 0.0;
        for(int i = 0; i < x.size(); i++) {
            C += fabs(TD[i]);
        }
        C = 1.0/C*elements.size();

        std::vector<std::vector<int> > nodetoglobal2 = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        int TDEGREE = Renumbering(nodetoglobal2);

        LILCSR<double> T = LILCSR<double>(TDEGREE, TDEGREE);
        std::vector<double> Y = std::vector<double>(TDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Te;
            Vector<double> Ye;
            LevelSet<double, ShapeFunction4Square, Gauss4Square>(Te, Ye, nodetoelement, elements[i], { 0 }, x, dt, tau, TD[i], phi);
            Assembling(T, Y, phi, Te, Ye, nodetoglobal2, nodetoelement, elements[i]);
		}

        CSR<double> Tmod = CSR<double>(T);	
        std::vector<double> result2 = ScalingCG(Tmod, Y, 100000, 1.0e-10);
        Disassembling(phi, result2, nodetoglobal2);

        
    }
    
    
    return 0;
}