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
	//******************************二次元熱伝導の熱伝導マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > HeatTransferTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _alpha, T _t) {
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Bマトリクス.....
		std::vector<std::vector<T> > B = std::vector<std::vector<T> >(2, std::vector<T>(3));
		B[0][0] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[0][1] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[0][2] = _nodes[_element[0]][1] - _nodes[_element[1]][1];
		B[1][0] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[1][1] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[1][2] = _nodes[_element[1]][0] - _nodes[_element[0]][0];
		B /= (2.0*A);

		return _alpha * Transpose(B)*B*_t;
	}


	//******************************二次元熱伝導の熱容量マトリクスを生成******************************
	template<class T>
	std::vector<std::vector<T> > HeatCapacityTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _rho, T _c, T _t) {
		//.....要素面積.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Ceマトリクス.....
		std::vector<std::vector<T> > Ce = std::vector<std::vector<T> >(3, std::vector<T>(3));
		Ce[0][0] = 1.0 / 6.0;	Ce[0][1] = 1.0 / 12.0;	Ce[0][2] = 1.0 / 12.0;
		Ce[1][0] = 1.0 / 12.0;	Ce[1][1] = 1.0 / 6.0;	Ce[1][2] = 1.0 / 12.0;
		Ce[2][0] = 1.0 / 12.0;	Ce[2][1] = 1.0 / 12.0;	Ce[2][2] = 1.0 / 6.0;
		Ce *= (_rho * _c*A*_t);
	
		return Ce;
	}
}