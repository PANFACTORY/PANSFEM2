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
	//*************************************************************************
	//	Direct method
	//*************************************************************************


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


	//******************************Get element stiffness matrix for Navier-Stokes equation******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _ubar, T _rho, T _mu) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> ubar = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			ubar(i, 0) = _ubar[_elementu[i]](0); ubar(i, 1) = _ubar[_elementu[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
		
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

			Vector<T> Ubar = ubar.Transpose()*M;		//	Advection velocity

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    K(i, j) = _rho*(M(i)*dMdX(0, j)*Ubar(0) + M(i)*dMdX(1, j)*Ubar(1)) + _mu*(2.0*dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));
					K(i, j + m) = _mu*dMdX(1, i)*dMdX(0, j);
					K(i + m, j) = _mu*dMdX(0, i)*dMdX(1, j);														
					K(i + m, j + m) = _rho*(M(i)*dMdX(1, j)*Ubar(1) + M(i)*dMdX(0, j)*Ubar(0)) + _mu*(dMdX(0, i)*dMdX(0, j) + 2.0*dMdX(1, i)*dMdX(1, j));	
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element consistent mass matrix for Navier-Stokes equation******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesConsistentMass(Matrix<T>& _Me, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _rho) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Me = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);

            Matrix<T> Me = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    Me(i, j) = _rho*M(i)*M(j);												
					Me(i + m, j + m) = _rho*M(i)*M(j);	
                }
            }
			_Me += Me*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element stiffness matrix for continuity expression******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void ContinuityStiffness(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ce = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();

			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

            Matrix<T> C = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
					C(i, j + 2*m) = -dMdX(0, i)*N(j);
					C(i + m, j + 2*m) = -dMdX(1, i)*N(j);
                    C(j + 2*m, i) = N(j)*dMdX(0, i);                                         						
					C(j + 2*m, i + m) = N(j)*dMdX(1, i);
                }
            }
			_Ce += C*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element stiffness matrix for Navier-Stokes equation with SUPG/PSPG******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesSUPGPSPGStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _ubar, T _rho, T _mu, T _dt) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> ubar = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			ubar(i, 0) = _ubar[_elementu[i]](0); ubar(i, 1) = _ubar[_elementu[i]](1);
		}

		T area = T();
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			area += J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
		
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

			Vector<T> Ubar = ubar.Transpose()*M;		//	Advection velocity

			T UNorm = Ubar.Norm();
			T he = T();
			if(UNorm > 1.0e-10) {
				Vector<T> dMdXU = dMdX.Transpose()*Ubar;
				for(int i = 0; i < m; i++) {
					he += fabs(dMdXU(i));
				}
				he = 2.0*UNorm/he;
			} else {
				he = sqrt(area);
			}
			T tau = 1.0/sqrt(pow(2.0/_dt, 2.0) + pow(2.0*UNorm/he, 2.0) + pow(4.0*_mu/(_rho*pow(he, 2.0)), 2.0));

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    K(i, j) = _rho*(Ubar(0)*Ubar(0)*dMdX(0, i)*dMdX(0, j) + Ubar(1)*Ubar(0)*dMdX(1, i)*dMdX(0, j) + Ubar(0)*Ubar(1)*dMdX(0, i)*dMdX(1, j) + Ubar(1)*Ubar(1)*dMdX(1, i)*dMdX(1, j));
					K(i + m, j + m) = _rho*(Ubar(0)*Ubar(0)*dMdX(0, i)*dMdX(0, j) + Ubar(1)*Ubar(0)*dMdX(1, i)*dMdX(0, j) + Ubar(0)*Ubar(1)*dMdX(0, i)*dMdX(1, j) + Ubar(1)*Ubar(1)*dMdX(1, i)*dMdX(1, j));	
                }
				for(int j = 0; j < n; j++) {
					K(j + 2*m, i) = _rho*(Ubar(0)*dNdX(0, j)*dMdX(0, i) + Ubar(1)*dNdX(0, j)*dMdX(1, i));
					K(j + 2*m, i + m) = _rho*(Ubar(0)*dNdX(1, j)*dMdX(0, i) + Ubar(1)*dNdX(1, j)*dMdX(1, i));
				}
            }
			_Ke += tau*K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element consistent mass matrix for Navier-Stokes equation with PSPG/SUPG******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesSUPGPSPGConsistentMass(Matrix<T>& _Me, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _ubar, T _rho, T _mu, T _dt) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Me = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> ubar = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			ubar(i, 0) = _ubar[_elementu[i]](0); ubar(i, 1) = _ubar[_elementu[i]](1);
		}

		T area = T();
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			area += J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

			Vector<T> Ubar = ubar.Transpose()*M;		//	Advection velocity

			T UNorm = Ubar.Norm();
			T he = T();
			if(UNorm > 1.0e-10) {
				Vector<T> dMdXU = dMdX.Transpose()*Ubar;
				for(int i = 0; i < m; i++) {
					he += fabs(dMdXU(i));
				}
				he = 2.0*UNorm/he;
			} else {
				he = sqrt(area);
			}
			T tau = 1.0/sqrt(pow(2.0/_dt, 2.0) + pow(2.0*UNorm/he, 2.0) + pow(4.0*_mu/(_rho*pow(he, 2.0)), 2.0));

            Matrix<T> Me = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    Me(i, j) = _rho*(Ubar(0)*dMdX(0, i)*M(j) + Ubar(1)*dMdX(1, i)*M(j));												
					Me(i + m, j + m) = _rho*(Ubar(0)*dMdX(0, i)*M(j) + Ubar(1)*dMdX(1, i)*M(j));
                }
				for(int j = 0; j < n; j++) {
					Me(j + 2*m, i) = _rho*dNdX(0, j)*M(i);
					Me(j + 2*m, i + m) = _rho*dNdX(1, j)*M(i);
				}
            }
			_Me += tau*Me*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element stiffness matrix for continuity expression with SUPG/PSPG******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void ContinuitySUPGPSPGStiffness(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _ubar, T _rho, T _mu, T _dt) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ce = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		Matrix<T> ubar = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			ubar(i, 0) = _ubar[_elementu[i]](0); ubar(i, 1) = _ubar[_elementu[i]](1);
		}

		T area = T();
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			area += J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
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

			Vector<T> Ubar = ubar.Transpose()*M;		//	Advection velocity

			T UNorm = Ubar.Norm();
			T he = T();
			if(UNorm > 1.0e-10) {
				Vector<T> dMdXU = dMdX.Transpose()*Ubar;
				for(int i = 0; i < m; i++) {
					he += fabs(dMdXU(i));
				}
				he = 2.0*UNorm/he;
			} else {
				he = sqrt(area);
			}
			T tau = 1.0/sqrt(pow(2.0/_dt, 2.0) + pow(2.0*UNorm/he, 2.0) + pow(4.0*_mu/(_rho*pow(he, 2.0)), 2.0));

            Matrix<T> C = Matrix<T>(2*m + n, 2*m + n);
			for(int i = 0; i < n; i++){
            	for(int j = 0; j < m; j++){
					C(j, i + 2*m) = Ubar(0)*dMdX(0, j)*dNdX(0, i) + Ubar(1)*dMdX(1, j)*dNdX(0, i);
					C(j + m, i + 2*m) = Ubar(0)*dMdX(0, j)*dNdX(1, i) + Ubar(1)*dMdX(1, j)*dNdX(1, i);
                }
				for(int j = 0; j < n; j++) {
					C(i + 2*m, j + 2*m) = dNdX(0, i)*dNdX(0, j) + dNdX(1, i)*dNdX(1, j);   
				}
            }
			_Ce += tau*C*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element body force vector******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC, class F>
	void NavierStokesBodyForce(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, F _f) {
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

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

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Vector<T> x = Xp.Transpose()*N;
			Vector<T> b = _f(x);

            Vector<T> Fe = Vector<T>(2*m + n);
            for(int i = 0; i < m; i++) {
                Fe(i) = b(0)*M(i);	Fe(i + m) = b(1)*M(i);
            }

			_Fe += Fe*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//*************************************************************************
	//	Decoupled method
	//*************************************************************************


	//******************************Get element consistent mass matrix for decoupled method******************************
	template<class T, template<class>class SF, template<class>class IC>
	void NavierStokesDecoupledConsistentMass(Matrix<T>& _Me, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _rho) {
		assert(_doulist.size() == 2);

		int m = _element.size();   //  Number of shapefunction for velosity u

		_Me = Matrix<T>(2*m, 2*m);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], m + i);
		}

		Matrix<T> X = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			X(i, 0) = _x[_element[i]](0); X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

            Matrix<T> M = Matrix<T>(2*m, 2*m);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    M(i, j) = _rho*N(i)*N(j);	M(i, j + m) = T();
                    M(i + m, j) = T();			M(i + m, j + m) = _rho*N(i)*N(j);
                }
            }
			_Me += M*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element lumped mass matrix for decoupled method******************************
	template<class T, template<class>class SF, template<class>class IC>
	void NavierStokesDecoupledLumpedMass(Matrix<T>& _Me, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _rho) {
		assert(_doulist.size() == 2);

		int m = _element.size();   //  Number of shapefunction for velosity u

		_Me = Matrix<T>(2*m, 2*m);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], m + i);
		}

		Matrix<T> X = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			X(i, 0) = _x[_element[i]](0); X(i, 1) = _x[_element[i]](1);
		}

		T Area = T();
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();            
			Area += J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}

		for(int i = 0; i < m; i++){
			_Me(i, i) = _rho*Area/(T)_element.size();
			_Me(i + m, i + m) = _rho*Area/(T)_element.size();
		}
	}


	//******************************Get element matrix and vector to get auxiliary velocity******************************
	template<class T, template<class>class SF, template<class>class IC>
	void NavierStokesAuxiliaryVelocity(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, T _rho, T _mu) {
		assert(_doulist.size() == 2);

		int m = _element.size();   //  Number of shapefunction for velosity u

		_Fe = Vector<T>(2*m);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], m + i);
		}

		Matrix<T> X = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			X(i, 0) = _x[_element[i]](0); X(i, 1) = _x[_element[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _u[_element[i]](0); u(i, 1) = _u[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> U = u.Transpose()*N;
			Matrix<T> dUdX = dNdX*u;

			Vector<T> F = Vector<T>(2*m);
            for(int i = 0; i < m; i++){
				F(i)     = _rho*N(i)*(U(0)*dUdX(0, 0) + U(1)*dUdX(1, 0)) + _mu*(dNdX(0, i)*dUdX(0, 0) + dNdX(1, i)*dUdX(1, 0));
				F(i + m) = _rho*N(i)*(U(0)*dUdX(0, 1) + U(1)*dUdX(1, 1)) + _mu*(dNdX(0, i)*dUdX(0, 1) + dNdX(1, i)*dUdX(1, 1));
            }
			_Fe -= F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element matrix and vector of Pressure-Poisson equation******************************
	template<class T, template<class>class SF, template<class>class IC>
	void NavierStokesPressurePoisson(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, T _rho, T _dt) {
		assert(_doulist.size() == 1);

		int m = _element.size();   //  Number of shapefunction for velosity u

		_Ke = Matrix<T>(m, m);
		_Fe = Vector<T>(m);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < m; i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			X(i, 0) = _x[_element[i]](0); X(i, 1) = _x[_element[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _u[_element[i]](0); u(i, 1) = _u[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> dUdX = dNdX*u;

            Matrix<T> K = Matrix<T>(m, m);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    K(i, j) = dNdX(0, i)*dNdX(0, j) + dNdX(1, i)*dNdX(1, j);
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];

			Vector<T> F = Vector<T>(m);
            for(int i = 0; i < m; i++){
				F(i) = -_rho/_dt*N(i)*(dUdX(0, 0) + dUdX(1, 1));
            }
			_Fe += F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element matrix and vector to get next step velocity******************************
	template<class T, template<class>class SF, template<class>class IC>
	void NavierStokesNextstepVelocity(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<Vector<T> >& _p, T _rho) {
		assert(_doulist.size() == 2);

		int m = _element.size();   //  Number of shapefunction for velosity u

		_Fe = Vector<T>(2*m);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], m + i);
		}

		Matrix<T> X = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			X(i, 0) = _x[_element[i]](0); X(i, 1) = _x[_element[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _u[_element[i]](0); u(i, 1) = _u[_element[i]](1);
		}

		Vector<T> p = Vector<T>(m);
		for(int i = 0; i < m; i++){
			p(i) = _p[_element[i]](0);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> U = u.Transpose()*N;
			Vector<T> dPdX = dNdX*p;

			Vector<T> F = Vector<T>(2*m);
            for(int i = 0; i < m; i++){
				F(i)     = N(i)*dPdX(0);
				F(i + m) = N(i)*dPdX(1);
            }
			_Fe -= F*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//*************************************************************************
	//	Adjoint equation
	//*************************************************************************


	//******************************Get element matrix of adjoint equation for NS******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesAdjoint(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementv, const std::vector<int>& _elementv, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementq, const std::vector<int>& _elementq, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _up, T _rho, T _mu){
		assert(_doulist.size() == 3);

		int m = _elementv.size();   //  Number of shapefunction for velosity u
        int n = _elementq.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementv = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementv[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementv[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementq = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementq[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

		Matrix<T> Xv = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xv(i, 0) = _x[_elementv[i]](0); Xv(i, 1) = _x[_elementv[i]](1);
		}

        Matrix<T> Xq = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xq(i, 0) = _x[_elementq[i]](0); Xq(i, 1) = _x[_elementq[i]](1);
		}

		Matrix<T> u = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			u(i, 0) = _up[_elementv[i]](0); u(i, 1) = _up[_elementv[i]](1);
		}

		Vector<T> p = Vector<T>(n);
		for(int i = 0; i < n; i++){
			p(i) = _up[_elementq[i]](2);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xq;
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
                    K(i, j) = _rho*(M(i)*M(j)*dUdX(0, 0) + dMdX(0, i)*U(0)*M(j) + dMdX(1, i)*U(1)*M(j)) + _mu*(2.0*dMdX(0, i)*dMdX(0, j) + dMdX(1, i)*dMdX(1, j));
					K(i, j + m) = _rho*M(i)*M(j)*dUdX(0, 1) + _mu*dMdX(1, i)*dMdX(0, j);
                    K(i + m, j) = _rho*M(i)*M(j)*dUdX(1, 0) + _mu*dMdX(0, i)*dMdX(1, j);
					K(i + m, j + m) = _rho*(dMdX(0, i)*U(0)*M(j) + M(i)*M(j)*dUdX(1, 1) + dMdX(1, i)*U(1)*M(j)) + _mu*(dMdX(0, i)*dMdX(0, j) + 2.0*dMdX(1, i)*dMdX(1, j));	
                }
                for(int j = 0; j < n; j++){
					K(i, j + 2*m) = -dMdX(0, i)*N(j);
					K(i + m, j + 2*m) = -dMdX(1, i)*N(j);
                    K(j + 2*m, i) = N(j)*dMdX(0, i);                                         						
					K(j + 2*m, i + m) = N(j)*dMdX(1, i);
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element body force matrix with Brinkman model******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void NavierStokesBrinkmanBodyForce(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, Matrix<T> _k) {
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    K(i, j) = M(i)*_k(0, 0)*M(j);		K(i, j + m) = M(i)*_k(0, 1)*M(j);
                    K(i + m, j) = M(i)*_k(1, 0)*M(j);   K(i + m, j + m) = M(i)*_k(1, 1)*M(j);  
                }
            }
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}