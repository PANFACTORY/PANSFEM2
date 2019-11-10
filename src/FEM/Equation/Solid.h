//*****************************************************************************
//Title		:src/FEM/Equation/Solid.h
//Author	:Tanabe Yuta
//Date		:2019/10/10
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../Controller/ShapeFunction.h"
#include "../Controller/IntegrationConstant.h"


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************Linear Isotropic Elastic Solid 3D********************
	template<class T>
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

		//----------Get Gauss integration constant----------
		std::vector<std::vector<T> > GP = GP3D8<T>;
		std::vector<std::vector<T> > GW = GW3D8<T>;
		
		//----------Loop of Gauss Integration----------
		for (int g = 0; g < 8; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = dNdr8I3<T>(GP[g]);

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
			_Ke = _Ke + B.Transpose()*C*B*J*GW[g][0] * GW[g][1] * GW[g][2];
		}
	}


	//********************Total Lagrange Solid 3D********************
	template<class T>
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

		//----------Get Gauss integration constant---------
		std::vector<std::vector<T> > GP = GP3D8<T>;
		std::vector<std::vector<T> > GW = GW3D8<T>;

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < 8; g++) {
			//----------Get difference of shape function by r----------
			Matrix<T> dNdr = dNdr8I3<T>(GP[g]);

			//----------Get difference of shape function by X----------
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse() * dNdr;

			//----------Get deformation gradient----------
			Matrix<T> Z = (dNdX * U).Transpose();

			//----------Generate initial displacement stiffness matrix----------
			Matrix<T> I = Matrix<T>(3, 3);
			I(0, 0)  = 1.0;	I(0, 1) = 0.0;	I(0, 2) = 0.0;
			I(1, 0)  = 0.0;	I(1, 1) = 1.0;	I(1, 2) = 0.0;
			I(2, 0)  = 0.0;	I(2, 1) = 0.0;	I(2, 2) = 1.0;
			Matrix<T> F = I + Z;
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
			_Ke += BL.Transpose()*C*BL*J*GW[g][0] * GW[g][1] * GW[g][2];

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

			Matrix<T> S = Matrix<T>(9, 9);
			S(0, 0) = Sv(0);	S(0, 1) = 0.0;		S(0, 2) = 0.0;			S(0, 3) = Sv(3);	S(0, 4) = 0.0;		S(0, 5) = 0.0;			S(0, 6) = Sv(5);	S(0, 7) = 0.0;		S(0, 8) = 0.0;
			S(1, 0) = 0.0;		S(1, 1) = Sv(0);	S(1, 2) = 0.0;			S(1, 3) = 0.0;		S(1, 4) = Sv(3);	S(1, 5) = 0.0;			S(1, 6) = 0.0;		S(1, 7) = Sv(5);	S(1, 8) = 0.0;
			S(2, 0) = 0.0;		S(2, 1) = 0.0;		S(2, 2) = Sv(0);		S(2, 3) = 0.0;		S(2, 4) = 0.0;		S(2, 5) = Sv(3);		S(2, 6) = 0.0;		S(2, 7) = 0.0;		S(2, 8) = Sv(5);

			S(3, 0) = Sv(3);	S(3, 1) = 0.0;		S(3, 2) = 0.0;			S(3, 3) = Sv(1);	S(3, 4) = 0.0;		S(3, 5) = 0.0;			S(3, 6) = Sv(4);	S(3, 7) = 0.0;		S(3, 8) = 0.0;
			S(4, 0) = 0.0;		S(4, 1) = Sv(3);	S(4, 2) = 0.0;			S(4, 3) = 0.0;		S(4, 4) = Sv(1);	S(4, 5) = 0.0;			S(4, 6) = 0.0;		S(4, 7) = Sv(4);	S(4, 8) = 0.0;
			S(5, 0) = 0.0;		S(5, 1) = 0.0;		S(5, 2) = Sv(3);		S(5, 3) = 0.0;		S(5, 4) = 0.0;		S(5, 5) = Sv(1);		S(5, 6) = 0.0;		S(5, 7) = 0.0;		S(5, 8) = Sv(4);

			S(6, 0) = Sv(5);	S(6, 1) = 0.0;		S(6, 2) = 0.0;			S(6, 3) = Sv(4);	S(6, 4) = 0.0;		S(6, 5) = 0.0;			S(6, 6) = Sv(2);	S(6, 7) = 0.0;		S(6, 8) = 0.0;
			S(7, 0) = 0.0;		S(7, 1) = Sv(5);	S(7, 2) = 0.0;			S(7, 3) = 0.0;		S(7, 4) = Sv(4);	S(7, 5) = 0.0;			S(7, 6) = 0.0;		S(7, 7) = Sv(2);	S(7, 8) = 0.0;
			S(8, 0) = 0.0;		S(8, 1) = 0.0;		S(8, 2) = Sv(5);		S(8, 3) = 0.0;		S(8, 4) = 0.0;		S(8, 5) = Sv(4);		S(8, 6) = 0.0;		S(8, 7) = 0.0;		S(8, 8) = Sv(2);

			_Ke += BNL.Transpose()*S*BNL*J*GW[g][0] * GW[g][1] * GW[g][2];

			//----------Update load vector----------
			_Fe += BL.Transpose()*Sv*J*GW[g][0] * GW[g][1] * GW[g][2];
		}
	}
}