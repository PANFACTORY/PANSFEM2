//*****************************************************************************
//Title		:src/FEM/Equation/Stokes.h
//Author	:Tanabe Yuta
//Date		:2020/02/05
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //******************************Get element advection matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void Stokes(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp, T _mu) {
		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p
        
        //----------Initialize element matrix----------
		_Ke = Matrix<T>(2*m + n, 2*m + n);

		//----------Generate cordinate matrix X of velosity u----------
		Matrix<T> Xu = Matrix<T>(0, 2);
		for(int i = 0; i < m; i++){
			Xu = Xu.Vstack(_x[_elementu[i]].Transpose());
		}

        //----------Generate cordinate matrix X of pressure p----------
		Matrix<T> Xp = Matrix<T>(0, 2);
		for(int i = 0; i < n; i++){
			Xp = Xp.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
            //----------Get shape function for pressure p----------
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Get shape function for velocity u----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

            //----------Get K matrix----------
            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    K(3*i, 3*j) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));  K(3*i, 3*j + 1) = T();                                                      K(3*i, 3*j + 2) = M(i)*dNdX(0, j);
                    K(3*i + 1, 3*j) = T();                                              K(3*i + 1, 3*j + 1) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));  K(3*i + 1, 3*j + 2) = M(i)*dNdX(1, j);
                    K(3*i + 2, 3*j) = M(j)*dNdX(0, i);                                  K(3*i + 2, 3*j + 1) = M(j)*dNdX(1, i);                                      K(3*i + 2, 3*j + 2) = T();
                }
                for(int j = n; j < m; j++){
                    K(3*i, n + 2*j) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));   K(3*i, n + 2*j + 1) = T();
                    K(3*i + 1, n + 2*j) = T();                                               K(3*i + 1, n + 2*j + 1) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));
                }
            }
            for(int i = n; i < m; i++){
                for(int j = 0; j < n; j++){
                    K(n + 2*i, 3*j) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));   K(n + 2*i, 3*j + 1) = T();
                    K(n + 2*i + 1, 3*j) = T();                                               K(n + 2*i + 1, 3*j + 1) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));
                }
                for(int j = n; j < m; j++){
                    K(n + 2*i, n + 2*j) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));   K(n + 2*i, n + 2*j + 1) = T();
                    K(n + 2*i + 1, n + 2*j) = T();                                               K(n + 2*i + 1, n + 2*j + 1) = _mu*(dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));
                }
            }
			
			//----------Update element matrix----------
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}