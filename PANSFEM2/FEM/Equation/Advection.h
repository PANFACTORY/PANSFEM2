//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/HeatTransfer.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/LAOperation.h"


namespace PANSFEM2 {
	//******************************二次元移流方程式の移流マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > AdvectionTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _cx, T _cy, T _t) {
		//.....移流マトリクスの確保.....
		std::vector<std::vector<T> > Se = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Bマトリクス.....
		std::vector<std::vector<T> > B = std::vector<std::vector<T> >(2, std::vector<T>(3));
		B[0][0] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[0][1] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[0][2] = _nodes[_element[0]][1] - _nodes[_element[1]][1];
		B[1][0] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[1][1] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[1][2] = _nodes[_element[1]][0] - _nodes[_element[0]][0];
		B = B / (2.0*A);

		//.....Seマトリクス.....
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Se[i][j] = (_cx*B[0][j] + _cy * B[1][j])*A / 3.0;
			}
		}
		
		return Se;
	}


	//******************************二次元質量マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > MassTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _t) {
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Meマトリクス.....
		std::vector<std::vector<T> > Me = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		Me[0][0] = 1.0 / 6.0;	Me[0][1] = 1.0 / 12.0;	Me[0][2] = 1.0 / 12.0;
		Me[1][0] = 1.0 / 12.0;	Me[1][1] = 1.0 / 6.0;	Me[1][2] = 1.0 / 12.0;
		Me[2][0] = 1.0 / 12.0;	Me[2][1] = 1.0 / 12.0;	Me[2][2] = 1.0 / 6.0;
		Me = Me * (A*_t);

		return Me;
	}
}