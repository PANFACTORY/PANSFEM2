//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/PlaneStrain.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Point.h"
#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
	//******************************平面ひずみモデル三角形一次要素******************************
	template<class T>
	void PlaneStrainTri(LILCSR<T>& _K, std::vector<T>& _F, std::vector<Point<T> >& _nodes, std::vector<int>& _nodestoelement, T _E, T _V, T _t) {
		//.....節点座標.....
		Point<double> p0 = _nodes[_nodestoelement[0]];
		Point<double> p1 = _nodes[_nodestoelement[1]];
		Point<double> p2 = _nodes[_nodestoelement[2]];
		
		//.....要素面積.....
		double A = 0.5*((p0.x[0] - p2.x[0])*(p1.x[1] - p2.x[1]) - (p2.x[1] - p0.x[1])*(p2.x[0] - p1.x[0]));

		//.....Bマトリクス.....
		double B[3][6];
		B[0][0] = 0.5*(p1.x[1] - p2.x[1]) / A;	B[0][1] = 0.0;							B[0][2] = 0.5*(p2.x[1] - p0.x[1]) / A;	B[0][3] = 0.0;							B[0][4] = 0.5*(p0.x[1] - p1.x[1]) / A;	B[0][5] = 0.0;
		B[1][0] = 0.0;							B[1][1] = 0.5*(p2.x[0] - p1.x[0]) / A;	B[1][2] = 0.0;							B[1][3] = 0.5*(p0.x[0] - p2.x[0]) / A;	B[1][4] = 0.0;							B[1][5] = 0.5*(p1.x[0] - p0.x[0]) / A;
		B[2][0] = 0.5*(p2.x[0] - p1.x[0]) / A;	B[2][1] = 0.5*(p1.x[1] - p2.x[1]) / A;	B[2][2] = 0.5*(p0.x[0] - p2.x[0]) / A;	B[2][3] = 0.5*(p2.x[1] - p0.x[1]) / A;	B[2][4] = 0.5*(p1.x[0] - p0.x[0]) / A;	B[2][5] = 0.5*(p0.x[1] - p1.x[1]) / A;

		//.....Dマトリクス.....
		double coef = _E / ((1.0 - 2.0*_V)*(1.0 + _V));
		double D[3][3];
		D[0][0] = coef * (1.0 - _V);	D[0][1] = coef * _V;			D[0][2] = 0.0;
		D[1][0] = D[0][1];				D[1][1] = coef * (1.0 - _V);	D[1][2] = 0.0;
		D[2][0] = D[0][2];				D[2][1] = D[1][2];				D[2][2] = coef * 0.5*(1.0 - 2.0*_V);

		//.....DBマトリクス.....
		double DB[3][6];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				DB[i][j] = 0.0;
				for (int k = 0; k < 3; k++) {
					DB[i][j] += D[i][k] * B[k][j];
				}
			}
		}

		//.....Ke(=∫BtDBdv)マトリクス.....
		double Ke[6][6];
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				Ke[i][j] = 0.0;
				for (int k = 0; k < 3; k++) {
					Ke[i][j] += B[k][i] * DB[k][j];
				}
				Ke[i][j] *= A * _t;
			}
		}

		//.....アセンブリング.....
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				_K.set(_nodestoelement[i] * 2, _nodestoelement[j] * 2, _K.get(_nodestoelement[i] * 2, _nodestoelement[j] * 2) + Ke[i * 2][j * 2]);
				_K.set(_nodestoelement[i] * 2, _nodestoelement[j] * 2 + 1, _K.get(_nodestoelement[i] * 2, _nodestoelement[j] * 2 + 1) + Ke[i * 2][j * 2 + 1]);
				_K.set(_nodestoelement[i] * 2 + 1, _nodestoelement[j] * 2, _K.get(_nodestoelement[i] * 2 + 1, _nodestoelement[j] * 2) + Ke[i * 2 + 1][j * 2]);
				_K.set(_nodestoelement[i] * 2 + 1, _nodestoelement[j] * 2 + 1, _K.get(_nodestoelement[i] * 2 + 1, _nodestoelement[j] * 2 + 1) + Ke[i * 2 + 1][j * 2 + 1]);
			}
		}
	}
}