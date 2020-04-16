//*****************************************************************************
//Title		:src/FEM/Equation/Solid.h
//Author	:Tanabe Yuta
//Date		:2019/10/10
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************Linear Isotropic Elastic Solid 3D********************
	template<class T, template<class>class SF, template<class>class IC>
	void LinearIsotropicElasticSolid(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _E, T _V) {
		//----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(3 * _element.size(), 3 * _element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(_element.size(), 3);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 3; j++){
				X(i, j) = _x[_element[i]](j);
			}
		}

		//----------Generate D matrix----------
		Matrix<T> C = Matrix<T>(6, 6);
		C(0, 0) = 1.0 - _V;	C(0, 1) = _V;		C(0, 2) = _V;		C(0, 3) = 0.0;					C(0, 4) = 0.0;					C(0, 5) = 0.0;
		C(1, 0) = _V;		C(1, 1) = 1.0 - _V;	C(1, 2) = _V;		C(1, 3) = 0.0;					C(1, 4) = 0.0;					C(1, 5) = 0.0;
		C(2, 0) = _V;		C(2, 1) = _V;		C(2, 2) = 1.0 - _V;	C(2, 3) = 0.0;					C(2, 4) = 0.0;					C(2, 5) = 0.0;
		C(3, 0) = 0.0;		C(3, 1) = 0.0;		C(3, 2) = 0.0;		C(3, 3) = 0.5*(1.0 - 2.0*_V);	C(3, 4) = 0.0;					C(3, 5) = 0.0;
		C(4, 0) = 0.0;		C(4, 1) = 0.0;		C(4, 2) = 0.0;		C(4, 3) = 0.0;					C(4, 4) = 0.5*(1.0 - 2.0*_V);	C(4, 5) = 0.0;
		C(5, 0) = 0.0;		C(5, 1) = 0.0;		C(5, 2) = 0.0;		C(5, 3) = 0.0;					C(5, 4) = 0.0;					C(5, 5) = 0.5*(1.0 - 2.0*_V);
		C *= _E / ((1.0 + _V)*(1.0 - 2.0*_V));

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse() * dNdr;

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(6, 3 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 3 * n) = dNdX(0, n);	B(0, 3 * n + 1) = 0.0;			B(0, 3 * n + 2) = 0.0;
				B(1, 3 * n) = 0.0;			B(1, 3 * n + 1) = dNdX(1, n);	B(1, 3 * n + 2) = 0.0;
				B(2, 3 * n) = 0.0;			B(2, 3 * n + 1) = 0.0;			B(2, 3 * n + 2) = dNdX(2, n);
				B(3, 3 * n) = dNdX(1, n);	B(3, 3 * n + 1) = dNdX(0, n);	B(3, 3 * n + 2) = 0.0;
				B(4, 3 * n) = 0.0;			B(4, 3 * n + 1) = dNdX(2, n);	B(4, 3 * n + 2) = dNdX(1, n);
				B(5, 3 * n) = dNdX(2, n);	B(5, 3 * n + 1) = 0.0;			B(5, 3 * n + 2) = dNdX(0, n);
			}

			//----------Update element stiffness matrix----------
			_Ke = _Ke + B.Transpose()*C*B*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1] * IC<T>::Weights[g][2];
		}
	}


	//********************Total Lagrange Solid 3D********************
	template<class T, template<class>class SF, template<class>class IC>
	void TotalLagrangeSolid(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<int>& _element, T _E, T _V) {
		//----------Initialize element stiffness matrix and load vector----------
		_Ke = Matrix<T>(3 * _element.size(), 3 * _element.size());
		_Fe = Vector<T>(3 * _element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(_element.size(), 3);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 3; j++){
				X(i, j) = _x[_element[i]](j);
			}
		}

		//----------Generate displacement matrix U----------
		Matrix<T> U = Matrix<T>(_element.size(), 3);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 3; j++){
				U(i, j) = _u[_element[i]](j);
			}
		}

		//----------Generate D matrix----------
		Matrix<T> C = Matrix<T>(6, 6);
		C(0, 0) = 1.0 - _V;	C(0, 1) = _V;		C(0, 2) = _V;		C(0, 3) = 0.0;					C(0, 4) = 0.0;					C(0, 5) = 0.0;
		C(1, 0) = _V;		C(1, 1) = 1.0 - _V;	C(1, 2) = _V;		C(1, 3) = 0.0;					C(1, 4) = 0.0;					C(1, 5) = 0.0;
		C(2, 0) = _V;		C(2, 1) = _V;		C(2, 2) = 1.0 - _V;	C(2, 3) = 0.0;					C(2, 4) = 0.0;					C(2, 5) = 0.0;
		C(3, 0) = 0.0;		C(3, 1) = 0.0;		C(3, 2) = 0.0;		C(3, 3) = 0.5*(1.0 - 2.0*_V);	C(3, 4) = 0.0;					C(3, 5) = 0.0;
		C(4, 0) = 0.0;		C(4, 1) = 0.0;		C(4, 2) = 0.0;		C(4, 3) = 0.0;					C(4, 4) = 0.5*(1.0 - 2.0*_V);	C(4, 5) = 0.0;
		C(5, 0) = 0.0;		C(5, 1) = 0.0;		C(5, 2) = 0.0;		C(5, 3) = 0.0;					C(5, 4) = 0.0;					C(5, 5) = 0.5*(1.0 - 2.0*_V);
		C *= _E / ((1.0 + _V)*(1.0 - 2.0*_V));

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function by r----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function by X----------
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse() * dNdr;

			//----------Get deformation gradient----------
			Matrix<T> Z = (dNdX * U).Transpose();

			//----------Generate initial displacement stiffness matrix----------
			Matrix<T> F = Identity<T>(3) + Z;
			Matrix<T> BL = Matrix<T>(6, 3 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				BL(0, 3 * n) = F(0, 0) * dNdX(0, n);						BL(0, 3 * n + 1) = F(1, 0) * dNdX(0, n);						BL(0, 3 * n + 2) = F(2, 0) * dNdX(0, n);
				BL(1, 3 * n) = F(0, 1) * dNdX(1, n);						BL(1, 3 * n + 1) = F(1, 1) * dNdX(1, n);						BL(1, 3 * n + 2) = F(2, 1) * dNdX(1, n);
				BL(2, 3 * n) = F(0, 2) * dNdX(2, n);						BL(2, 3 * n + 1) = F(1, 2) * dNdX(2, n);						BL(2, 3 * n + 2) = F(2, 2) * dNdX(2, n);
				BL(3, 3 * n) = F(0, 1) * dNdX(0, n) + F(0, 0) * dNdX(1, n);	BL(3, 3 * n + 1) = F(1, 1) * dNdX(0, n) + F(1, 0) * dNdX(1, n);	BL(3, 3 * n + 2) = F(2, 1) * dNdX(0, n) + F(2, 0) * dNdX(1, n);
				BL(4, 3 * n) = F(0, 2) * dNdX(1, n) + F(0, 1) * dNdX(2, n);	BL(4, 3 * n + 1) = F(1, 2) * dNdX(1, n) + F(1, 1) * dNdX(2, n);	BL(4, 3 * n + 2) = F(2, 2) * dNdX(1, n) + F(2, 1) * dNdX(2, n);
				BL(5, 3 * n) = F(0, 0) * dNdX(2, n) + F(0, 2) * dNdX(0, n);	BL(5, 3 * n + 1) = F(1, 0) * dNdX(2, n) + F(1, 2) * dNdX(0, n);	BL(5, 3 * n + 2) = F(2, 0) * dNdX(2, n) + F(2, 2) * dNdX(0, n);
			}

			//----------Update element stiffness matrix----------
			_Ke += BL.Transpose()*C*BL*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1] * IC<T>::Weights[g][2];

			//----------Get Green-Lagrange strain----------
			Matrix<T> E = (Z + Z.Transpose() + Z.Transpose()*Z) / 2.0;

			//----------Get Piola-Kirchhoff stress----------
			Vector<T> Ev = Vector<T>({ E(0, 0), E(1, 1), E(2, 2), E(0, 1) + E(1, 0), E(1, 2) + E(2, 1), E(2, 0) + E(0, 2) });
			Vector<T> Sv = C * Ev;

			//----------Generate initial stress stiffness matrix----------
			Matrix<T> BNL = Matrix<T>(9, 3 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				BNL(0, 3 * n) = dNdX(0, n);	BNL(0, 3 * n + 1) = 0.0;		BNL(0, 3 * n + 2) = 0.0;
				BNL(1, 3 * n) = 0.0;		BNL(1, 3 * n + 1) = dNdX(0, n);	BNL(1, 3 * n + 2) = 0.0;
				BNL(2, 3 * n) = 0.0;		BNL(2, 3 * n + 1) = 0.0;		BNL(2, 3 * n + 2) = dNdX(0, n);
				BNL(3, 3 * n) = dNdX(1, n);	BNL(3, 3 * n + 1) = 0.0;		BNL(3, 3 * n + 2) = 0.0;
				BNL(4, 3 * n) = 0.0;		BNL(4, 3 * n + 1) = dNdX(1, n);	BNL(4, 3 * n + 2) = 0.0;
				BNL(5, 3 * n) = 0.0;		BNL(5, 3 * n + 1) = 0.0;		BNL(5, 3 * n + 2) = dNdX(1, n);
				BNL(6, 3 * n) = dNdX(2, n);	BNL(6, 3 * n + 1) = 0.0;		BNL(6, 3 * n + 2) = 0.0;
				BNL(7, 3 * n) = 0.0;		BNL(7, 3 * n + 1) = dNdX(2, n);	BNL(7, 3 * n + 2) = 0.0;
				BNL(8, 3 * n) = 0.0;		BNL(8, 3 * n + 1) = 0.0;		BNL(8, 3 * n + 2) = dNdX(2, n);
			}

			Matrix<T> S00 = Sv(0)*Identity<T>(3);	Matrix<T> S01 = Sv(3)*Identity<T>(3);	Matrix<T> S02 = Sv(5)*Identity<T>(3);
			Matrix<T> S10 = Sv(3)*Identity<T>(3);	Matrix<T> S11 = Sv(1)*Identity<T>(3);	Matrix<T> S12 = Sv(4)*Identity<T>(3);
			Matrix<T> S20 = Sv(5)*Identity<T>(3);	Matrix<T> S21 = Sv(4)*Identity<T>(3);	Matrix<T> S22 = Sv(2)*Identity<T>(3);
			Matrix<T> S = (S00.Hstack(S01.Hstack(S02))).Vstack((S10.Hstack(S11.Hstack(S12))).Vstack((S20.Hstack(S21.Hstack(S22)))));

			_Ke += BNL.Transpose()*S*BNL*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1] * IC<T>::Weights[g][2];

			//----------Update load vector----------
			_Fe += BL.Transpose()*Sv*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1] * IC<T>::Weights[g][2];
		}
	}


	//********************Updated Lagrange Solid 3D********************
	template<class T, template<class>class SF, template<class>class IC>
	void UpdatedLagrangeSolid(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<int>& _element, T _E, T _V) {
		//----------Initialize element stiffness matrix and load vector----------
		_Ke = Matrix<T>(3*_element.size(), 3*_element.size());
		_Fe = Vector<T>(3*_element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 3);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Generate displacement matrix U----------
		Matrix<T> U = Matrix<T>(0, 3);
		for(int i = 0; i < _element.size(); i++){
			U = U.Vstack(_u[_element[i]].Transpose());
		}

		//----------Generate cordinate matrix x----------
		Matrix<T> x = X + U;

		//----------Get Lame parameters----------
		T mu0 = 0.5*_E/(1.0 + _V);
		T lambda0 = 2.0*mu0*_V/(1.0 - 2.0*_V);

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function by r----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dxdr = dNdr*x;
			T J = dxdr.Determinant();
			Matrix<T> dNdx = dxdr.Inverse()*dNdr;
			Matrix<T> F = ((dNdr*X).Inverse()*dNdr*x).Transpose();	

			//----------Generate initial displacement stiffness matrix----------
			T detF = F.Determinant();
			T mu = (mu0 - lambda0*log(detF))/detF;
			T lambda = lambda0/detF;
			Matrix<T> C = Matrix<T>(6, 6);
			C(0, 0) = lambda + 2.0*mu;	C(0, 1) = lambda;			C(0, 2) = lambda;			C(0, 3) = T();	C(0, 4) = T();	C(0, 5) = T();
			C(1, 0) = lambda;			C(1, 1) = lambda + 2.0*mu;	C(1, 2) = lambda;			C(1, 3) = T();	C(1, 4) = T();	C(1, 5) = T();
			C(2, 0) = lambda;			C(2, 1) = lambda;			C(2, 2) = lambda + 2.0*mu;	C(2, 3) = T();	C(2, 4) = T();	C(2, 5) = T();
			C(3, 0) = T();				C(3, 1) = T();				C(3, 2) = T();				C(3, 3) = mu;	C(3, 4) = T();	C(3, 5) = T();
			C(4, 0) = T();				C(4, 1) = T();				C(4, 2) = T();				C(4, 3) = T();	C(4, 4) = mu;	C(4, 5) = T();
			C(5, 0) = T();				C(5, 1) = T();				C(5, 2) = T();				C(5, 3) = T();	C(5, 4) = T();	C(5, 5) = mu;
			Matrix<T> BL = Matrix<T>(6, 3*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				BL(0, 3*n) = dNdx(0, n);	BL(0, 3*n + 1) = T();			BL(0, 3*n + 2) = T();
				BL(1, 3*n) = T();			BL(1, 3*n + 1) = dNdx(1, n);	BL(1, 3*n + 2) = T();
				BL(2, 3*n) = T();			BL(2, 3*n + 1) = T();			BL(2, 3*n + 2) = dNdx(2, n);
				BL(3, 3*n) = dNdx(1, n);	BL(3, 3*n + 1) = dNdx(0, n);	BL(3, 3*n + 2) = T();
				BL(4, 3*n) = dNdx(2, n);	BL(4, 3*n + 1) = T();			BL(4, 3*n + 2) = dNdx(0, n);
				BL(5, 3*n) = T();			BL(5, 3*n + 1) = dNdx(2, n);	BL(5, 3*n + 2) = dNdx(1, n);
			}
			_Ke += BL.Transpose()*C*BL*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1]*IC<T>::Weights[g][2];

			//----------Generate initial stress stiffness matrix----------
			Matrix<T> S = (mu0/detF)*(F*F.Transpose() - Identity<T>(3)) + (lambda0*log(detF)/detF)*Identity<T>(3);
			Matrix<T> BNL = Matrix<T>(9, 3*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				BNL(0, 3*n) = dNdx(0, n);	BNL(0, 3*n + 1) = T();			BNL(0, 3*n + 2) = T();
				BNL(1, 3*n) = T();			BNL(1, 3*n + 1) = dNdx(0, n);	BNL(1, 3*n + 2) = T();
				BNL(2, 3*n) = T();			BNL(2, 3*n + 1) = T();			BNL(2, 3*n + 2) = dNdx(0, n);
				BNL(3, 3*n) = dNdx(1, n);	BNL(3, 3*n + 1) = T();			BNL(3, 3*n + 2) = T();
				BNL(4, 3*n) = T();			BNL(4, 3*n + 1) = dNdx(1, n);	BNL(4, 3*n + 2) = T();
				BNL(5, 3*n) = T();			BNL(5, 3*n + 1) = T();			BNL(5, 3*n + 2) = dNdx(1, n);
				BNL(6, 3*n) = dNdx(2, n);	BNL(6, 3*n + 1) = T();			BNL(6, 3*n + 2) = T();
				BNL(7, 3*n) = T();			BNL(7, 3*n + 1) = dNdx(2, n);	BNL(7, 3*n + 2) = T();
				BNL(8, 3*n) = T();			BNL(8, 3*n + 1) = T();			BNL(8, 3*n + 2) = dNdx(2, n);
			}
			Matrix<T> S00 = S(0, 0)*Identity<T>(3);	Matrix<T> S01 = S(0, 1)*Identity<T>(3);	Matrix<T> S02 = S(0, 2)*Identity<T>(3);
			Matrix<T> S10 = S(1, 0)*Identity<T>(3);	Matrix<T> S11 = S(1, 1)*Identity<T>(3);	Matrix<T> S12 = S(1, 2)*Identity<T>(3);
			Matrix<T> S20 = S(2, 0)*Identity<T>(3);	Matrix<T> S21 = S(2, 1)*Identity<T>(3);	Matrix<T> S22 = S(2, 2)*Identity<T>(3);
			Matrix<T> Sm = (S00.Hstack(S01.Hstack(S02))).Vstack((S10.Hstack(S11.Hstack(S12))).Vstack((S20.Hstack(S21.Hstack(S22)))));
			_Ke += BNL.Transpose()*Sm*BNL*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1]*IC<T>::Weights[g][2];

			//----------Update load vector----------
			Vector<T> Sv = Vector<T>({ S(0, 0), S(1, 1), S(2, 2), S(0, 1), S(0, 2), S(1, 2) });	
			_Fe += BL.Transpose()*Sv*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1]*IC<T>::Weights[g][2];
		}
	}










	//********************Linear Isotropic Elastic Solid 3D********************
	template<class T, template<class>class SF, template<class>class IC>
	void SolidLinearIsotropicElastic(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _E, T _V) {
		assert(_doulist.size() == 3);

		_Ke = Matrix<T>(3*_element.size(), 3*_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(3));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 3*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 3*i + 1);
			_nodetoelement[i][2] = std::make_pair(_doulist[2], 3*i + 2);
		}
		
		Matrix<T> X = Matrix<T>(_element.size(), 3);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);	X(i, 2) = _x[_element[i]](2);
		}

		Matrix<T> C = Matrix<T>(6, 6);
		C(0, 0) = 1.0 - _V;	C(0, 1) = _V;		C(0, 2) = _V;		C(0, 3) = T();					C(0, 4) = T();					C(0, 5) = T();
		C(1, 0) = _V;		C(1, 1) = 1.0 - _V;	C(1, 2) = _V;		C(1, 3) = T();					C(1, 4) = T();					C(1, 5) = T();
		C(2, 0) = _V;		C(2, 1) = _V;		C(2, 2) = 1.0 - _V;	C(2, 3) = T();					C(2, 4) = T();					C(2, 5) = T();
		C(3, 0) = T();		C(3, 1) = T();		C(3, 2) = T();		C(3, 3) = 0.5*(1.0 - 2.0*_V);	C(3, 4) = T();					C(3, 5) = T();
		C(4, 0) = T();		C(4, 1) = T();		C(4, 2) = T();		C(4, 3) = T();					C(4, 4) = 0.5*(1.0 - 2.0*_V);	C(4, 5) = T();
		C(5, 0) = T();		C(5, 1) = T();		C(5, 2) = T();		C(5, 3) = T();					C(5, 4) = T();					C(5, 5) = 0.5*(1.0 - 2.0*_V);
		C *= _E/((1.0 + _V)*(1.0 - 2.0*_V));

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> B = Matrix<T>(6, 3 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 3*n) = dNdX(0, n);	B(0, 3*n + 1) = T();		B(0, 3*n + 2) = T();
				B(1, 3*n) = T();		B(1, 3*n + 1) = dNdX(1, n);	B(1, 3*n + 2) = T();
				B(2, 3*n) = T();		B(2, 3*n + 1) = T();		B(2, 3*n + 2) = dNdX(2, n);
				B(3, 3*n) = dNdX(1, n);	B(3, 3*n + 1) = dNdX(0, n);	B(3, 3*n + 2) = T();
				B(4, 3*n) = T();		B(4, 3*n + 1) = dNdX(2, n);	B(4, 3*n + 2) = dNdX(1, n);
				B(5, 3*n) = dNdX(2, n);	B(5, 3*n + 1) = T();		B(5, 3*n + 2) = dNdX(0, n);
			}

			_Ke += B.Transpose()*C*B*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1]*IC<T>::Weights[g][2];
		}
	}	
}