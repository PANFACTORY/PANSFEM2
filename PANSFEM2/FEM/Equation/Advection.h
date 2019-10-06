//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/HeatTransfer.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************二次元移流方程式の移流マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > AdvectionTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _cx, T _cy, T _t) {
		//.....移流マトリクスの確保.....
		std::vector<std::vector<T> > Se = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		
		//.....要素面積.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);
		
		//.....Bマトリクス.....
		T B[2][3];
		B[0][0] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[0][1] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[0][2] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;
		B[1][0] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[1][1] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[1][2] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;

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
	std::vector<std::vector<T> > MassTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _t) {
		//.....質量マトリクスの確保.....
		std::vector<std::vector<T> > Me = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));

		//.....要素面積.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);

		//.....Meマトリクス.....
		Me[0][0] = A*_t / 6.0;	Me[0][1] = A*_t / 12.0;	Me[0][2] = A*_t / 12.0;
		Me[1][0] = A*_t / 12.0;	Me[1][1] = A*_t / 6.0;	Me[1][2] = A*_t / 12.0;
		Me[2][0] = A*_t / 12.0;	Me[2][1] = A*_t / 12.0;	Me[2][2] = A*_t / 6.0;

		return Me;
	}
}