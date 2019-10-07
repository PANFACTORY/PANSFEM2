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
	//******************************�񎟌��ڗ��������̈ڗ��}�g���N�X�𐶐�******************************
	template<class T>
	std::vector<std::vector<T> > AdvectionTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _cx, T _cy, T _t) {
		//.....�ڗ��}�g���N�X�̊m��.....
		std::vector<std::vector<T> > Se = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		
		//.....�v�f�ʐ�.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....B�}�g���N�X.....
		std::vector<std::vector<T> > B = std::vector<std::vector<T> >(2, std::vector<T>(3));
		B[0][0] = _nodes[_element[1]][1] - _nodes[_element[2]][1];	B[0][1] = _nodes[_element[2]][1] - _nodes[_element[0]][1];	B[0][2] = _nodes[_element[0]][1] - _nodes[_element[1]][1];
		B[1][0] = _nodes[_element[2]][0] - _nodes[_element[1]][0];	B[1][1] = _nodes[_element[0]][0] - _nodes[_element[2]][0];	B[1][2] = _nodes[_element[1]][0] - _nodes[_element[0]][0];
		B = B / (2.0*A);

		//.....Se�}�g���N�X.....
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Se[i][j] = (_cx*B[0][j] + _cy * B[1][j])*A / 3.0;
			}
		}
		
		return Se;
	}


	//******************************�񎟌����ʃ}�g���N�X�𐶐�******************************
	template<class T>
	std::vector<std::vector<T> > MassTri(std::vector<std::vector<T> >& _nodes, std::vector<int>& _element, T _t) {
		//.....�v�f�ʐ�.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Me�}�g���N�X.....
		std::vector<std::vector<T> > Me = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		Me[0][0] = 1.0 / 6.0;	Me[0][1] = 1.0 / 12.0;	Me[0][2] = 1.0 / 12.0;
		Me[1][0] = 1.0 / 12.0;	Me[1][1] = 1.0 / 6.0;	Me[1][2] = 1.0 / 12.0;
		Me[2][0] = 1.0 / 12.0;	Me[2][1] = 1.0 / 12.0;	Me[2][2] = 1.0 / 6.0;
		Me = Me * (A*_t);

		return Me;
	}
}