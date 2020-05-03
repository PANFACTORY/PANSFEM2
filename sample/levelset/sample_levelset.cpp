#include <iostream>
#include <vector>


#include "../../src/PrePost/Mesher/SquareMesh.h"


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
    

    //----------最適化ループ----------
    for(int t = 0; t < 200; t++) {
        //----------ステップ2：固定設計領域の有限要素離散化と数値解析----------

        //----------ステップ3：目的汎関数と制約汎関数の計算----------

        //----------ステップ4：収束条件の判定----------

        //----------ステップ5：目的汎関数と制約汎関数の設計感度の計算----------

        //----------ステップ6：レベルセット関数の更新，ステップ2に戻る----------
    }
    
    
    return 0;
}