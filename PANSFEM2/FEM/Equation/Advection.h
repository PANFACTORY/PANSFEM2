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
	//******************************�񎟌��ڗ��������̈ڗ��}�g���N�X�𐶐�******************************
	template<class T>
	std::vector<std::vector<T> > AdvectionTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _cx, T _cy, T _t) {
		//.....�ڗ��}�g���N�X�̊m��.....
		std::vector<std::vector<T> > Se = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));
		
		//.....�v�f�ʐ�.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);
		
		//.....B�}�g���N�X.....
		T B[2][3];
		B[0][0] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[0][1] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[0][2] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;
		B[1][0] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[1][1] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[1][2] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;

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
	std::vector<std::vector<T> > MassTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _t) {
		//.....���ʃ}�g���N�X�̊m��.....
		std::vector<std::vector<T> > Me = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));

		//.....�v�f�ʐ�.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);

		//.....Me�}�g���N�X.....
		Me[0][0] = A*_t / 6.0;	Me[0][1] = A*_t / 12.0;	Me[0][2] = A*_t / 12.0;
		Me[1][0] = A*_t / 12.0;	Me[1][1] = A*_t / 6.0;	Me[1][2] = A*_t / 12.0;
		Me[2][0] = A*_t / 12.0;	Me[2][1] = A*_t / 12.0;	Me[2][2] = A*_t / 6.0;

		return Me;
	}
}