//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/Solid.h
//Author	:Tanabe Yuta
//Date		:2019/10/10
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LAOperation.h"
#include "../Controller/ShapeFunction.h"
#include "../Controller/IntegrationConstant.h"


namespace PANSFEM2 {
	//********************Linear Isotropic Elastic Solid 3D********************
	template<class T>
	void LinearIsotropicElasticSolid3(std::vector<std::vector<T> >& _Ke, std::vector<std::vector<T> >& _x, std::vector<int>& _element, T _E, T _V) {
		//----------接線剛性行列，等価節点内力ベクトルの確保----------
		_Ke = std::vector<std::vector<T> >(3 * _element.size(), std::vector<T>(3 * _element.size(), T()));

		//----------座標行列を求める----------
		std::vector<std::vector<T> > X;
		for (auto nodeid : _element) {
			X.push_back(_x[nodeid]);
		}

		//----------構成則を求める----------
		std::vector<std::vector<T> > C = std::vector<std::vector<T> >(6, std::vector<T>(6));
		C[0][0] = 1.0 - _V;	C[0][1] = _V;		C[0][2] = _V;		C[0][3] = 0.0;					C[0][4] = 0.0;					C[0][5] = 0.0;
		C[1][0] = _V;		C[1][1] = 1.0 - _V;	C[1][2] = _V;		C[1][3] = 0.0;					C[1][4] = 0.0;					C[1][5] = 0.0;
		C[2][0] = _V;		C[2][1] = _V;		C[2][2] = 1.0 - _V;	C[2][3] = 0.0;					C[2][4] = 0.0;					C[2][5] = 0.0;
		C[3][0] = 0.0;		C[3][1] = 0.0;		C[3][2] = 0.0;		C[3][3] = 0.5*(1.0 - 2.0*_V);	C[3][4] = 0.0;					C[3][5] = 0.0;
		C[4][0] = 0.0;		C[4][1] = 0.0;		C[4][2] = 0.0;		C[4][3] = 0.0;					C[4][4] = 0.5*(1.0 - 2.0*_V);	C[4][5] = 0.0;
		C[5][0] = 0.0;		C[5][1] = 0.0;		C[5][2] = 0.0;		C[5][3] = 0.0;					C[5][4] = 0.0;					C[5][5] = 0.5*(1.0 - 2.0*_V);
		C *= _E / ((1.0 + _V)*(1.0 - 2.0*_V));

		//----------ガウス積分用の定数----------
		std::vector<std::vector<T> > GP = GP3D8<T>;
		std::vector<std::vector<T> > GW = GW3D8<T>;
		
		//----------数値積分ループ内----------
		for (int g = 0; g < 8; g++) {
			//----------形状関数行列を求める----------
			std::vector<std::vector<T> > dNdr = dNdr20I3<T>(GP[g]);

			//----------座標関連の計算----------
			std::vector<std::vector<T> > dXdr = dNdr * X;
			T J = Determinant3(dXdr);
			std::vector<std::vector<T> > dNdX = Inverse3(dXdr) * dNdr;

			//----------変位勾配テンソルFからひずみ―変位関係行列Bを求める----------
			std::vector<std::vector<T> > B = std::vector<std::vector<T> >(6, std::vector<T>(3 * _element.size()));
			for (int n = 0; n < _element.size(); n++) {
				B[0][3 * n] = dNdX[0][n];	B[0][3 * n + 1] = 0.0;			B[0][3 * n + 2] = 0.0;
				B[1][3 * n] = 0.0;			B[1][3 * n + 1] = dNdX[1][n];	B[1][3 * n + 2] = 0.0;
				B[2][3 * n] = 0.0;			B[2][3 * n + 1] = 0.0;			B[2][3 * n + 2] = dNdX[2][n];
				B[3][3 * n] = dNdX[1][n];	B[3][3 * n + 1] = dNdX[0][n];	B[3][3 * n + 2] = 0.0;
				B[4][3 * n] = 0.0;			B[4][3 * n + 1] = dNdX[2][n];	B[4][3 * n + 2] = dNdX[1][n];
				B[5][3 * n] = dNdX[2][n];	B[5][3 * n + 1] = 0.0;			B[5][3 * n + 2] = dNdX[0][n];
			}

			//----------初期変位マトリクスを求める----------
			_Ke += Transpose(B)*C*B*J*GW[g][0] * GW[g][1] * GW[g][2];
		}
	}


	//********************Total Lagrange Solid 3D********************
	template<class T>
	void TotalLagrangeSolid3(std::vector<std::vector<T> >& _Ke, std::vector<T>& _Fe, std::vector<std::vector<T> >& _x, std::vector<std::vector<T> >& _u, std::vector<int>& _element, T _E, T _V) {
		//----------接線剛性行列，等価節点内力ベクトルの確保----------
		_Ke = std::vector<std::vector<T> >(3 * _element.size(), std::vector<T>(3 * _element.size(), T()));
		_Fe = std::vector<T>(3 * _element.size(), T());

		//----------座標行列を求める----------
		std::vector<std::vector<T> > X;
		for (auto nodeid : _element) {
			X.push_back(_x[nodeid]);
		}

		//----------変位行列を求める----------
		std::vector<std::vector<T> > U;
		for (auto nodeid : _element) {
			U.push_back(_u[nodeid]);
		}

		//----------構成則を求める----------
		std::vector<std::vector<T> > C = std::vector<std::vector<T> >(6, std::vector<T>(6));
		C[0][0] = 1.0 - _V;	C[0][1] = _V;		C[0][2] = _V;		C[0][3] = 0.0;					C[0][4] = 0.0;					C[0][5] = 0.0;
		C[1][0] = _V;		C[1][1] = 1.0 - _V;	C[1][2] = _V;		C[1][3] = 0.0;					C[1][4] = 0.0;					C[1][5] = 0.0;
		C[2][0] = _V;		C[2][1] = _V;		C[2][2] = 1.0 - _V;	C[2][3] = 0.0;					C[2][4] = 0.0;					C[2][5] = 0.0;
		C[3][0] = 0.0;		C[3][1] = 0.0;		C[3][2] = 0.0;		C[3][3] = 0.5*(1.0 - 2.0*_V);	C[3][4] = 0.0;					C[3][5] = 0.0;
		C[4][0] = 0.0;		C[4][1] = 0.0;		C[4][2] = 0.0;		C[4][3] = 0.0;					C[4][4] = 0.5*(1.0 - 2.0*_V);	C[4][5] = 0.0;
		C[5][0] = 0.0;		C[5][1] = 0.0;		C[5][2] = 0.0;		C[5][3] = 0.0;					C[5][4] = 0.0;					C[5][5] = 0.5*(1.0 - 2.0*_V);
		C *= _E / ((1.0 + _V)*(1.0 - 2.0*_V));

		//----------ガウス積分用の定数---------
		std::vector<std::vector<T> > GP = GP3D8<T>;
		std::vector<std::vector<T> > GW = GW3D8<T>;

		//----------数値積分ループ内----------
		for (int g = 0; g < 8; g++) {
			//----------形状関数行列を求める----------
			std::vector<std::vector<T> > dNdr = dNdr20I3<T>(GP[g]);

			//----------座標関連の計算----------
			std::vector<std::vector<T> > dXdr = dNdr * X;
			T J = Determinant3(dXdr);
			std::vector<std::vector<T> > dNdX = Inverse3(dXdr) * dNdr;

			//----------変位勾配テンソルZを求める----------
			std::vector<std::vector<T> > Z = Transpose(dNdX * U);

			//----------変位勾配テンソルZからGreen-LagrangeひずみEを求める----------
			std::vector<std::vector<T> > E = 0.5 * (Z + Transpose(Z) + Transpose(Z)*Z);

			//----------Green-LagrangeひずみEから第二Piola-Kirchhoff応力Sを求める----------
			std::vector<T> Ev = { E[0][0], E[1][1], E[2][2], E[0][1] + E[1][0], E[1][2] + E[2][1], E[2][0] + E[0][2] };
			std::vector<T> Sv = C * Ev;

			//----------変位勾配テンソルFからひずみ―変位関係行列Bを求める----------
			std::vector<std::vector<T> > I = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
			std::vector<std::vector<T> > F = I + Z;
			std::vector<std::vector<T> > BL = std::vector<std::vector<T> >(6, std::vector<T>(3 * _element.size()));
			for (int n = 0; n < _element.size(); n++) {
				BL[0][3 * n] = F[0][0] * dNdX[0][n];						BL[0][3 * n + 1] = F[1][0] * dNdX[0][n];						BL[0][3 * n + 2] = F[2][0] * dNdX[0][n];
				BL[1][3 * n] = F[0][1] * dNdX[1][n];						BL[1][3 * n + 1] = F[1][1] * dNdX[1][n];						BL[1][3 * n + 2] = F[2][1] * dNdX[1][n];
				BL[2][3 * n] = F[0][2] * dNdX[2][n];						BL[2][3 * n + 1] = F[1][2] * dNdX[2][n];						BL[2][3 * n + 2] = F[2][2] * dNdX[2][n];
				BL[3][3 * n] = F[0][1] * dNdX[0][n] + F[0][0] * dNdX[1][n];	BL[3][3 * n + 1] = F[1][1] * dNdX[0][n] + F[1][0] * dNdX[1][n];	BL[3][3 * n + 2] = F[2][1] * dNdX[0][n] + F[2][0] * dNdX[1][n];
				BL[4][3 * n] = F[0][2] * dNdX[1][n] + F[0][1] * dNdX[2][n];	BL[4][3 * n + 1] = F[1][2] * dNdX[1][n] + F[1][1] * dNdX[2][n];	BL[4][3 * n + 2] = F[2][2] * dNdX[1][n] + F[2][1] * dNdX[2][n];
				BL[5][3 * n] = F[0][0] * dNdX[2][n] + F[0][2] * dNdX[0][n];	BL[5][3 * n + 1] = F[1][0] * dNdX[2][n] + F[1][2] * dNdX[0][n];	BL[5][3 * n + 2] = F[2][0] * dNdX[2][n] + F[2][2] * dNdX[0][n];
			}

			//----------初期変位マトリクスを求める----------
			_Ke += Transpose(BL)*C*BL*J*GW[g][0] * GW[g][1] * GW[g][2];

			//----------初期応力マトリクスを求める----------
			std::vector<std::vector<T> > BNL = std::vector<std::vector<T> >(9, std::vector<T>(3 * _element.size()));
			for (int n = 0; n < _element.size(); n++) {
				BNL[0][3 * n] = dNdX[0][n];	BNL[0][3 * n + 1] = 0.0;		BNL[0][3 * n + 2] = 0.0;
				BNL[1][3 * n] = 0.0;		BNL[1][3 * n + 1] = dNdX[0][n];	BNL[1][3 * n + 2] = 0.0;
				BNL[2][3 * n] = 0.0;		BNL[2][3 * n + 1] = 0.0;		BNL[2][3 * n + 2] = dNdX[0][n];
				BNL[3][3 * n] = dNdX[1][n];	BNL[3][3 * n + 1] = 0.0;		BNL[3][3 * n + 2] = 0.0;
				BNL[4][3 * n] = 0.0;		BNL[4][3 * n + 1] = dNdX[1][n];	BNL[4][3 * n + 2] = 0.0;
				BNL[5][3 * n] = 0.0;		BNL[5][3 * n + 1] = 0.0;		BNL[5][3 * n + 2] = dNdX[1][n];
				BNL[6][3 * n] = dNdX[2][n];	BNL[6][3 * n + 1] = 0.0;		BNL[6][3 * n + 2] = 0.0;
				BNL[7][3 * n] = 0.0;		BNL[7][3 * n + 1] = dNdX[2][n];	BNL[7][3 * n + 2] = 0.0;
				BNL[8][3 * n] = 0.0;		BNL[8][3 * n + 1] = 0.0;		BNL[8][3 * n + 2] = dNdX[2][n];
			}

			std::vector<std::vector<T> > S = std::vector<std::vector<T> >(9, std::vector<T>(9));
			S[0][0] = Sv[0];	S[0][1] = 0.0;		S[0][2] = 0.0;			S[0][3] = Sv[3];	S[0][4] = 0.0;		S[0][5] = 0.0;			S[0][6] = Sv[5];	S[0][7] = 0.0;		S[0][8] = 0.0;
			S[1][0] = 0.0;		S[1][1] = Sv[0];	S[1][2] = 0.0;			S[1][3] = 0.0;		S[1][4] = Sv[3];	S[1][5] = 0.0;			S[1][6] = 0.0;		S[1][7] = Sv[5];	S[1][8] = 0.0;
			S[2][0] = 0.0;		S[2][1] = 0.0;		S[2][2] = Sv[0];		S[2][3] = 0.0;		S[2][4] = 0.0;		S[2][5] = Sv[3];		S[2][6] = 0.0;		S[2][7] = 0.0;		S[2][8] = Sv[5];

			S[3][0] = Sv[3];	S[3][1] = 0.0;		S[3][2] = 0.0;			S[3][3] = Sv[1];	S[3][4] = 0.0;		S[3][5] = 0.0;			S[3][6] = Sv[4];	S[3][7] = 0.0;		S[3][8] = 0.0;
			S[4][0] = 0.0;		S[4][1] = Sv[3];	S[4][2] = 0.0;			S[4][3] = 0.0;		S[4][4] = Sv[1];	S[4][5] = 0.0;			S[4][6] = 0.0;		S[4][7] = Sv[4];	S[4][8] = 0.0;
			S[5][0] = 0.0;		S[5][1] = 0.0;		S[5][2] = Sv[3];		S[5][3] = 0.0;		S[5][4] = 0.0;		S[5][5] = Sv[1];		S[5][6] = 0.0;		S[5][7] = 0.0;		S[5][8] = Sv[4];

			S[6][0] = Sv[5];	S[6][1] = 0.0;		S[6][2] = 0.0;			S[6][3] = Sv[4];	S[6][4] = 0.0;		S[6][5] = 0.0;			S[6][6] = Sv[2];	S[6][7] = 0.0;		S[6][8] = 0.0;
			S[7][0] = 0.0;		S[7][1] = Sv[5];	S[7][2] = 0.0;			S[7][3] = 0.0;		S[7][4] = Sv[4];	S[7][5] = 0.0;			S[7][6] = 0.0;		S[7][7] = Sv[2];	S[7][8] = 0.0;
			S[8][0] = 0.0;		S[8][1] = 0.0;		S[8][2] = Sv[5];		S[8][3] = 0.0;		S[8][4] = 0.0;		S[8][5] = Sv[4];		S[8][6] = 0.0;		S[8][7] = 0.0;		S[8][8] = Sv[2];

			_Ke += Transpose(BNL)*S*BNL*J*GW[g][0] * GW[g][1] * GW[g][2];

			//----------第二Piola-Kirchhoff応力Sとひずみ―変位関係行列Bから要素節点等価内力ベクトルを求める----------
			_Fe += Transpose(BL)*Sv*J*GW[g][0] * GW[g][1] * GW[g][2];
		}
	}


	//********************Updated Lagrange Solid 3D********************
	//template<class T>

}