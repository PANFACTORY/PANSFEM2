//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/PlaneStrain.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************平面ひずみモデルの剛性マトリクスを生成******************************
	template<class T>
	void PlaneStrainTri(std::vector<std::vector<T> >& _Ke, std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _E, T _V, T _t) {
		//.....要素剛性マトリクスの確保.....
		_Ke = std::vector<std::vector<double> >(6, std::vector<double>(6, T()));
		
		//.....要素面積.....
		double A = 0.5*((_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0])*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) - (_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1])*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]));

		//.....Bマトリクス.....
		double B[3][6];
		B[0][0] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[0][1] = 0.0;																B[0][2] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[0][3] = 0.0;																B[0][4] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;	B[0][5] = 0.0;
		B[1][0] = 0.0;																B[1][1] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[1][2] = 0.0;																B[1][3] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[1][4] = 0.0;																B[1][5] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;
		B[2][0] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[2][1] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[2][2] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[2][3] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[2][4] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;	B[2][5] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;

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
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 3; k++) {
					_Ke[i][j] += B[k][i] * DB[k][j];
				}
				_Ke[i][j] *= A * _t;
			}
		}
	}
}