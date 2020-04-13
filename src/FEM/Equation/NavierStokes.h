//*****************************************************************************
//  Title       :src/FEM/Equation/NavierStokes.h
//  Author      :Tanabe Yuta
//  Date        :2020/02/10
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //******************************Get element stiffness matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokes(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<T>& _p, std::vector<int>& _elementu, std::vector<int>& _elementp, T _nu) {
		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p
        
        //----------Initialize element matrix----------
		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_Fe = Vector<T>(2*m + n);

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

		//----------Generate velocity matrix U----------
		Matrix<T> u = Matrix<T>(0, 2);
		for(int i = 0; i < m; i++){
			u = u.Vstack(_u[_elementu[i]].Transpose());
		}

		//----------Generate pressure vector p----------
		Vector<T> p = Vector<T>(n);
		for(int i = 0; i < n; i++){
			p(i) = _p[_elementp[i]];
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

			//----------Get velocity u and derivative of u----------
			Vector<T> U = u.Transpose()*M;	//	ui
			Matrix<T> dUdX = dMdX*u;		//	dujdXi

			//----------Get pressure p and derivative of p----------
			T P = p*N;						//	p
			Vector<T> dpdX = dNdX*p;		//	dpdXi

			//----------Get diffusion term----------
			Matrix<T> Sxx = Matrix<T>(m, m);
			Matrix<T> Sxy = Matrix<T>(m, m);
			Matrix<T> Syx = Matrix<T>(m, m);
			Matrix<T> Syy = Matrix<T>(m, m);
			for(int i = 0; i < m; i++) {
				for(int j = 0; j < m; j++) {
					Sxx(i, j) = 2.0*_nu*dMdX(0, i)*dMdX(0, j) + _nu*dMdX(1, i)*dMdX(1, j);
					Sxy(i, j) = _nu*dMdX(1, i)*dMdX(0, j);
					Syx(i, j) = _nu*dMdX(0, i)*dMdX(1, j);
					Syy(i, j) = _nu*dMdX(0, i)*dMdX(0, j) + 2.0*_nu*dMdX(1, i)*dMdX(1, j);
				}
			}

			//----------Get convective term----------
			Matrix<T> Kx0 = Matrix<T>(m, m);
			Matrix<T> Kx1 = Matrix<T>(m, m);
			Matrix<T> Ky0 = Matrix<T>(m, m);
			Matrix<T> Ky1 = Matrix<T>(m, m);
			for(int i = 0; i < m; i++) {
				for(int j = 0; j < m; j++) {
					Kx0(i, j) = M(i)*M(j);
					Kx1(i, j) = M(i)*dMdX(0, j);
					Ky0(i, j) = M(i)*M(j);
					Ky1(i, j) = M(i)*dMdX(1, j);
				}
			}

            //----------Get Ke matrix----------
            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    K(3*i, 3*j) = Kx0(i, j)*dUdX(0, 0) + Kx1(i, j)*U(0) + Ky1(i, j)*U(1) + Sxx(i, j) + Sxy(i, j);	K(3*i, 3*j + 1) = Ky0(i, j)*dUdX(1, 0);																	K(3*i, 3*j + 2) = -N(j)*dMdX(0, i);
                    K(3*i + 1, 3*j) = Kx0(i, j)*dUdX(0, 1);															K(3*i + 1, 3*j + 1) = Ky0(i, j)*dUdX(1, 1) + Ky1(i, j)*U(1) + Kx1(i, j)*U(0) + Syx(i, j) + Syy(i, j);	K(3*i + 1, 3*j + 2) = -N(j)*dMdX(1, i);
                    K(3*i + 2, 3*j) = N(i)*dMdX(0, j);																K(3*i + 2, 3*j + 1) = N(i)*dMdX(1, j);																	K(3*i + 2, 3*j + 2) = T();
                }
                for(int j = n; j < m; j++){
                    K(3*i, n + 2*j) = Kx0(i, j)*dUdX(0, 0) + Kx1(i, j)*U(0) + Ky1(i, j)*U(1) + Sxx(i, j) + Sxy(i, j);	K(3*i, n + 2*j + 1) = Ky0(i, j)*dUdX(1, 0);
                    K(3*i + 1, n + 2*j) = Kx0(i, j)*dUdX(0, 1);															K(3*i + 1, n + 2*j + 1) = Ky0(i, j)*dUdX(1, 1) + Ky1(i, j)*U(1) + Kx1(i, j)*U(0) + Syx(i, j) + Syy(i, j);
                    K(3*i + 2, n + 2*j) = N(i)*dMdX(0, j);                                         						K(3*i + 2, n + 2*j + 1) = N(i)*dMdX(1, j);
                }
            }
            for(int i = n; i < m; i++){
                for(int j = 0; j < n; j++){
                    K(n + 2*i, 3*j) = Kx0(i, j)*dUdX(0, 0) + Kx1(i, j)*U(0) + Ky1(i, j)*U(1) + Sxx(i, j) + Sxy(i, j);	K(n + 2*i, 3*j + 1) = Ky0(i, j)*dUdX(1, 0);                                     							K(n + 2*i, 3*j + 2) = -N(j)*dMdX(0, i);
                    K(n + 2*i + 1, 3*j) = Kx0(i, j)*dUdX(0, 1);															K(n + 2*i + 1, 3*j + 1) = Ky0(i, j)*dUdX(1, 1) + Ky1(i, j)*U(1) + Kx1(i, j)*U(0) + Syx(i, j) + Syy(i, j);	K(n + 2*i + 1, 3*j + 2) = -N(j)*dMdX(1, i);
                }
                for(int j = n; j < m; j++){
                    K(n + 2*i, n + 2*j) = Kx0(i, j)*dUdX(0, 0) + Kx1(i, j)*U(0) + Ky1(i, j)*U(1) + Sxx(i, j) + Sxy(i, j);	K(n + 2*i, n + 2*j + 1) = Ky0(i, j)*dUdX(1, 0);
                    K(n + 2*i + 1, n + 2*j) = Kx0(i, j)*dUdX(0, 1);															K(n + 2*i + 1, n + 2*j + 1) = Ky0(i, j)*dUdX(1, 1) + Ky1(i, j)*U(1) + Kx1(i, j)*U(0) + Syx(i, j) + Syy(i, j);
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];

			//----------Get Fe vector----------
			Vector<T> F = Vector<T>(2*m + n);
            for(int i = 0; i < n; i++){
				F(3*i)     = M(i)*U(0)*dUdX(0, 0) + M(i)*U(1)*dUdX(1, 0) + 2.0*_nu*dMdX(0, i)*dUdX(0, 0) + _nu*dMdX(1, i)*dUdX(1, 0) + _nu*dMdX(1, i)*dUdX(0, 1) - dMdX(0, i)*P;
				F(3*i + 1) = M(i)*U(0)*dUdX(0, 1) + M(i)*U(1)*dUdX(1, 1) + _nu*dMdX(0, i)*dUdX(1, 0) + _nu*dMdX(0, i)*dUdX(0, 1) + 2.0*_nu*dMdX(1, i)*dUdX(1, 1) - dMdX(1, i)*P;
				F(3*i + 2) = N(i)*dUdX(0, 0) + N(i)*dUdX(1, 1);
            }
            for(int i = n; i < m; i++){
                F(n + 2*i)     = M(i)*U(0)*dUdX(0, 0) + M(i)*U(1)*dUdX(1, 0) + 2.0*_nu*dMdX(0, i)*dUdX(0, 0) + _nu*dMdX(1, i)*dUdX(1, 0) + _nu*dMdX(1, i)*dUdX(0, 1) - dMdX(0, i)*P;
				F(n + 2*i + 1) = M(i)*U(0)*dUdX(0, 1) + M(i)*U(1)*dUdX(1, 1) + _nu*dMdX(0, i)*dUdX(1, 0) + _nu*dMdX(0, i)*dUdX(0, 1) + 2.0*_nu*dMdX(1, i)*dUdX(1, 1) - dMdX(1, i)*P;
            }
			_Fe -= F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}