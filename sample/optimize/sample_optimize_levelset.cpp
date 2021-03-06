#include <iostream>
#include <vector>
#include <numeric>


#include "../../src/LinearAlgebra/Models/Vector.h"
#include "../../src/PrePost/Mesher/SquareMesh.h"
#include "../../src/FEM/Equation/PlaneStress.h"
#include "../../src/FEM/Equation/ReactionDiffusion.h"
#include "../../src/FEM/Equation/General.h"
#include "../../src/FEM/Controller/ShapeFunction.h"
#include "../../src/FEM/Controller/GaussIntegration.h"
#include "../../src/FEM/Controller/BoundaryCondition.h"
#include "../../src/FEM/Controller/Assembling.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/LinearAlgebra/Solvers/CG.h"
#include "../../src/PrePost/Export/ExportToVTK.h"


using namespace PANSFEM2;


int main() {
    //----------パラメータの設定----------
    double Vmax = 0.5;
    double tau = 2.0e-4;
    double E0 = 1.0;
    double Emin = 1.0e-4;
    double nu = 0.3;
    double nvol = 100;
    double dt = 0.1;
    double d = -0.02;
    double p = 4.0;
    int tmax = 200;

    double A1 = -1.5*(1.0 - nu)*(1.0 - 14.0*nu + 15.0*pow(nu, 2.0))*E0/((1.0 + nu)*(7.0 - 5.0*nu)*pow(1.0 - 2.0*nu, 2.0));
    double A2 = 7.5*(1.0 - nu)*E0/((1.0 + nu)*(7.0 - 5.0*nu));
    double c = A1/(A1 + 2.0*A2);


    //----------固定設計領域と境界条件の設定----------
    SquareMesh<double> mesh = SquareMesh<double>(60.0, 40.0, 60, 40);
    std::vector<Vector<double> > x = mesh.GenerateNodes();
    std::vector<std::vector<int> > elements = mesh.GenerateElements();
    std::vector<std::pair<std::pair<int, int>, double> > ufixed = mesh.GenerateFixedlist({ 0, 1 }, [](Vector<double> _x){
        if(abs(_x(0)) < 1.0e-5) {
            return true;
        }
        return false;
    });
    std::vector<std::pair<std::pair<int, int>, double> > qfixed = mesh.GenerateFixedlist({ 1 }, [](Vector<double> _x){
        if(abs(_x(0) - 60.0) < 1.0e-5 && abs(_x(1) - 20.0) < 1.0 +  1.0e-5) {
            return true;
        }
        return false;
    });
    for(auto& qfixedi : qfixed) {
        qfixedi.second = -1.0;
    }
    std::vector<std::pair<std::pair<int, int>, double> > phifixed = mesh.GenerateFixedlist({ 0 }, [](Vector<double> _x){
        if(abs(_x(0) - 0.0) < 1.0e-5 || abs(_x(0) - 60.0) < 1.0e-5 || abs(_x(1) - 0.0) < 1.0e-5 || abs(_x(1) - 40.0) < 1.0e-5) {
            return true;
        }
        return false;
    });


    std::vector<Vector<double> > phi = std::vector<Vector<double> >(x.size(), { 1.0 });     //  φ
    std::vector<double> str = std::vector<double>(elements.size(), 1.0);                    //  χ(φ)
    double volInit = std::accumulate(str.begin(), str.end(), 0.0)/(double)elements.size();
    std::vector<double> objective = std::vector<double>(tmax);
    

    //----------最適化ループ----------
    for(int t = 0; t < tmax; t++) {
        //----------固定設計領域の有限要素法解析----------
        std::vector<Vector<double> > u = std::vector<Vector<double> >(x.size(), Vector<double>(2));
        std::vector<std::vector<int> > nodetoglobal = std::vector<std::vector<int> >(x.size(), std::vector<int>(2, 0));
        
        SetDirichlet(u, nodetoglobal, ufixed);
        int KDEGREE = Renumbering(nodetoglobal);

        LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
        std::vector<double> F = std::vector<double>(KDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStressStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, Emin + str[i]*(E0 - Emin), nu, 1.0);
            Assembling(K, F, u, Ke, nodetoglobal, nodetoelement, elements[i]);
		}

        Assembling(F, qfixed, nodetoglobal);

        CSR<double> Kmod = CSR<double>(K);	
        std::vector<double> result = ScalingCG(Kmod, F, 100000, 1.0e-10);
        Disassembling(u, result, nodetoglobal);     


        //----------ファイル出力----------
        std::ofstream fout("sample/optimize/result" + std::to_string(t) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(x, fout);
		AddElementToVTK(elements, fout);
		AddElementTypes(std::vector<int>(elements.size(), 9), fout);
		AddPointVectors(u, "u", fout, true);
        AddPointScalers(phi, "phi", fout, false);
		AddElementScalers(str, "str", fout, true);
		fout.close();


        //----------トポロジカルデリバティブと体積の計算----------
        std::vector<Vector<double> > TD = std::vector<Vector<double> >(elements.size(), Vector<double>(1));
        std::vector<int> NC = std::vector<int>(x.size(), 0);
        for (int i = 0; i < elements.size(); i++) {
			std::vector<std::vector<std::pair<int, int> > > nodetoelement;
            Matrix<double> Ke;
            PlaneStressStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, Emin + str[i]*(E0 - Emin), nu, 1.0);
            Vector<double> ue = ElementVector(u, nodetoelement, elements[i]);
            objective[t] += ue*(Ke*ue);
            PlaneStressStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0, 1 }, x, (A1 + 2.0*A2)*(1.0 - pow(c, 2.0)), c, 1.0);
            TD[i](0) = (1.0e-4 + str[i]*(1.0 - 1.0e-4))*ue*(Ke*ue);
        }
        std::vector<Vector<double> > TDN = InterpolateNodalFromElemental<double, Vector>(x.size(), Vector<double>(1), TD, elements);
        
        double vol = std::accumulate(str.begin(), str.end(), 0.0)/(double)elements.size();
        double ex = Vmax + (volInit - Vmax)*std::max(0.0, 1.0 - (t + 1)/(double)nvol);
        double lambda = std::accumulate(TDN.begin(), TDN.end(), Vector<double>(1))(0)/(double)x.size()*exp(p*((vol - ex)/ex + d));


        //----------収束条件の判定----------
        if(t > nvol && fabs(vol - Vmax) < 0.005 && 
            fabs(objective[t] - objective[t - 5]) < 0.01*fabs(objective[t]) && 
            fabs(objective[t] - objective[t - 4]) < 0.01*fabs(objective[t]) && 
            fabs(objective[t] - objective[t - 3]) < 0.01*fabs(objective[t]) &&
            fabs(objective[t] - objective[t - 2]) < 0.01*fabs(objective[t]) &&
            fabs(objective[t] - objective[t - 1]) < 0.01*fabs(objective[t])) {
            
            std::cout << "----------Convergence----------" << std::endl;
            break;
        }

        std::cout << "t = " << t << "\tCompliance = " << objective[t]/(double)elements.size() << "\tVolume = " << vol << "\tLambda = " << lambda << std::endl;

        
        //----------レベルセット関数の更新----------
        double C = 0.0;
        for(int i = 0; i < x.size(); i++) {
            C += fabs(TDN[i](0));
        }
        C = elements.size()/C;

        std::vector<std::vector<int> > nodetoglobal2 = std::vector<std::vector<int> >(x.size(), std::vector<int>(1, 0));
        SetDirichlet(phi, nodetoglobal2, phifixed);
        int TDEGREE = Renumbering(nodetoglobal2);

        LILCSR<double> T = LILCSR<double>(TDEGREE, TDEGREE);
        std::vector<double> Y = std::vector<double>(TDEGREE, 0.0);

		for (int i = 0; i < elements.size(); i++) {
            std::vector<std::vector<std::pair<int, int> > > nodetoelement;
			Matrix<double> Me, Ke;
            Vector<double> Fe;
            ReactionDiffusionConsistentMass<double, ShapeFunction4Square, Gauss4Square>(Me, nodetoelement, elements[i], { 0 }, x);
            ReactionDiffusionStiffness<double, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0 }, x, tau*elements.size());
            ReactionDiffusionReaction<double, ShapeFunction4Square, Gauss4Square>(Fe, nodetoelement, elements[i], { 0 }, x, TDN, [&](double _u, Vector<double> _dudX) {
                return C*(_u - lambda);
            });
			Vector<double> phie = ElementVector(phi, nodetoelement, elements[i]);
			Matrix<double> Te = Me/dt + Ke;
			Vector<double> Ye = Me/dt*phie + Fe; 
			Assembling(T, Y, phi, Te, nodetoglobal2, nodetoelement, elements[i]);
            Assembling(Y, Ye, nodetoglobal2, nodetoelement, elements[i]);
		}

        CSR<double> Tmod = CSR<double>(T);	
        std::vector<double> result2 = ScalingCG(Tmod, Y, 100000, 1.0e-10);
        Disassembling(phi, result2, nodetoglobal2);

        for(int i = 0; i < x.size(); i++) {
            phi[i](0) = std::max(std::min(1.0, phi[i](0)), -1.0);
        }

        std::vector<Vector<double> > phie = InterpolateElementalFromNodal<double, Vector>(Vector<double>(1), phi, elements);
        for(int i = 0; i < elements.size(); i++) {
            if(phie[i](0) < 0.0) {
                str[i] = 0.0;
            } else {
                str[i] = 1.0;
            }
        }
    }
    
    
    return 0;
}