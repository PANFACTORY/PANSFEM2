//*****************************************************************************
//  Title       :	src/FEM/Equation/NavierStokes.h
//  Author      :	Tanabe Yuta
//  Date        :	2020/02/10
//  Copyright   :	(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Get element tengent matrix for Navier-Stokes equation******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesTangent(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _up, T _rho, T _mu) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_Fe = Vector<T>(2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

		Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _up[_elementu[i]](0); u(i, 1) = _up[_elementu[i]](1);
		}

		Vector<T> p = Vector<T>(n);
		for(int i = 0; i < n; i++){
			p(i) = _up[_elementp[i]](2);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

			Vector<T> U = u.Transpose()*M;
			Matrix<T> dUdX = dMdX*u;

			T P = p*N;
			Vector<T> dpdX = dNdX*p;

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    K(i, j) = _rho*M(i)*M(j)*dUdX(0, 0) + _rho*M(i)*dMdX(0, j)*U(0) + _rho*M(i)*dMdX(1, j)*U(1) + 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j) + _mu*dMdX(1, i)*dMdX(0, j);	
					K(i, j + m) = _rho*M(i)*M(j)*dUdX(1, 0);																	
                    K(i + m, j) = _rho*M(i)*M(j)*dUdX(0, 1);																		
					K(i + m, j + m) = _rho*M(i)*M(j)*dUdX(1, 1) + _rho*M(i)*dMdX(1, j)*U(1) + _rho*M(i)*dMdX(0, j)*U(0) + _mu*dMdX(0, i)*dMdX(1, j) + _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);	
                }
                for(int j = 0; j < n; j++){
					K(i, j + 2*m) = -dMdX(0, i)*N(j);
					K(i + m, j + 2*m) = -dMdX(1, i)*N(j);
                    K(j + 2*m, i) = N(j)*dMdX(0, i);                                         						
					K(j + 2*m, i + m) = N(j)*dMdX(1, i);
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];

			Vector<T> F = Vector<T>(2*m + n);
            for(int i = 0; i < m; i++){
				F(i)     = _rho*M(i)*U(0)*dUdX(0, 0) + _rho*M(i)*U(1)*dUdX(1, 0) + 2.0*_mu*dMdX(0, i)*dUdX(0, 0) + _mu*dMdX(1, i)*dUdX(1, 0) + _mu*dMdX(1, i)*dUdX(0, 1) - dMdX(0, i)*P;
				F(i + m) = _rho*M(i)*U(0)*dUdX(0, 1) + _rho*M(i)*U(1)*dUdX(1, 1) + _mu*dMdX(0, i)*dUdX(1, 0) + _mu*dMdX(0, i)*dUdX(0, 1) + 2.0*_mu*dMdX(1, i)*dUdX(1, 1) - dMdX(1, i)*P;
            }
            for(int i = 0; i < n; i++){
                F(i + 2*m) = N(i)*dUdX(0, 0) + N(i)*dUdX(1, 1);
            }
			_Fe -= F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element matrix and vector for dynamic Navier-Stokes equation******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesDynamic(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _up, T _rho, T _mu, T dt) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_Fe = Vector<T>(2*m + n); 
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

		Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _up[_elementu[i]](0); u(i, 1) = _up[_elementu[i]](1);
		}

		Vector<T> p = Vector<T>(n);
		for(int i = 0; i < n; i++){
			p(i) = _up[_elementp[i]](2);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

			Vector<T> U = u.Transpose()*M;
			Matrix<T> dUdX = dMdX*u;

			T P = p*N;
			Vector<T> dpdX = dNdX*p;

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    K(i, j) = _rho*M(i)*M(j)/dt + _rho*M(i)*dMdX(0, j)*U(0) + _rho*M(i)*dMdX(1, j)*U(1) + 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j) + _mu*dMdX(1, i)*dMdX(0, j);															
					K(i + m, j + m) = _rho*M(i)*M(j)/dt + _rho*M(i)*dMdX(1, j)*U(1) + _rho*M(i)*dMdX(0, j)*U(0) + _mu*dMdX(0, i)*dMdX(1, j) + _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);	
                }
                for(int j = 0; j < n; j++){
					K(i, j + 2*m) = -dMdX(0, i)*N(j);
					K(i + m, j + 2*m) = -dMdX(1, i)*N(j);
                    K(j + 2*m, i) = N(j)*dMdX(0, i);                                         						
					K(j + 2*m, i + m) = N(j)*dMdX(1, i);
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];

			Vector<T> F = Vector<T>(2*m + n);
			for(int i = 0; i < m; i++) {
				for(int j = 0; j < m; j++) {
					F(i) += _rho*M(i)*M(j)/dt*_up[_elementu[j]](0);
					F(i + m) += _rho*M(i)*M(j)/dt*_up[_elementu[j]](1);
				}
			}
			_Fe += F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}