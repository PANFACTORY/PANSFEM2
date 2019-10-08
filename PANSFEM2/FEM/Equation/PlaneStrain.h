//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/PlaneStrain.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LAOperation.h"


namespace PANSFEM2 {
	//******************************平面ひずみモデルの剛性マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > PlaneStrainTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _E, T _V, T _t) {
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Bマトリクス.....
		std::vector<std::vector<T> > B = std::vector<std::vector<T> >(3, std::vector<T>(6));
		B[0][0] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[0][1] = 0.0;												B[0][2] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[0][3] = 0.0;												B[0][4] = _nodes[_element[0]][1] - _nodes[_element[1]][1];	B[0][5] = 0.0;
		B[1][0] = 0.0;												B[1][1] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[1][2] = 0.0;												B[1][3] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[1][4] = 0.0;												B[1][5] = _nodes[_element[1]][0] - _nodes[_element[0]][0];
		B[2][0] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[2][1] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[2][2] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[2][3] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[2][4] = _nodes[_element[1]][0] - _nodes[_element[0]][0];	B[2][5] = _nodes[_element[0]][1] - _nodes[_element[1]][1];
		B /= (2.0*A);

		//.....Dマトリクス.....
		std::vector<std::vector<T> > D = std::vector<std::vector<T> >(3, std::vector<T>(3));
		D[0][0] = 1.0 - _V;	D[0][1] = _V;		D[0][2] = 0.0;
		D[1][0] = D[0][1];	D[1][1] = 1.0 - _V;	D[1][2] = 0.0;
		D[2][0] = D[0][2];	D[2][1] = D[1][2];	D[2][2] = 0.5*(1.0 - 2.0*_V);
		D *= (_E / ((1.0 - 2.0*_V)*(1.0 + _V)));

		//.....Keマトリクスを計算して返す.....
		return Transpose(B)*D*B*A*_t;
	}
}