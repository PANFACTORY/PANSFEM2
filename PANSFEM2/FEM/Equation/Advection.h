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
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Bマトリクス.....
		std::vector<std::vector<T> > B = std::vector<std::vector<T> >(2, std::vector<T>(3));
		B[0][0] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[0][1] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[0][2] = _nodes[_element[0]][1] - _nodes[_element[1]][1];
		B[1][0] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[1][1] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[1][2] = _nodes[_element[1]][0] - _nodes[_element[0]][0];
		B /= (2.0*A);

		//.....Nマトリクス.....
		std::vector<std::vector<T> > N = std::vector<std::vector<T> >(3, std::vector<T>(1));
		N[0][0] = A / 3.0;
		N[1][0] = A / 3.0;
		N[2][0] = A / 3.0;

		//.....移流速度マトリクス.....
		std::vector<T> c = { _cx, _cy };
		std::vector<std::vector<T> > C = Transpose(c);	

		return N*C*B*_t;
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
		Me *= (A*_t);

		return Me;
	}
}