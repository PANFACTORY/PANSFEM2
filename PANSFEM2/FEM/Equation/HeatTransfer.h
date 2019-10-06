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
	//******************************�񎟌��M�`���̔M�`���}�g���N�X�𐶐�******************************
	template<class T>
	std::vector<std::vector<T> > HeatTransferTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _alpha, T _t) {
		//.....�v�f�M�`���}�g���N�X�̊m��.....
		std::vector<std::vector<T> > Ke = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));

		//.....�v�f�ʐ�.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);

		//.....B�}�g���N�X.....
		T B[2][3];
		B[0][0] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[0][1] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[0][2] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;
		B[1][0] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[1][1] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[1][2] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;

		//.....Ke(=�烿BtBdv)�}�g���N�X.....
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 2; k++) {
					Ke[i][j] += B[k][i] * B[k][j];
				}
				Ke[i][j] *= _alpha * _t;
			}
		}

		return Ke;
	}


	//******************************�񎟌��M�`���̔M�e�ʃ}�g���N�X�𐶐�******************************
	template<class T>
	std::vector<std::vector<T> > HeatCapacityTri(std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _rho, T _c, T _t) {
		//.....�v�f�M�e�ʃ}�g���N�X�̊m��.....
		std::vector<std::vector<T> > Ce = std::vector<std::vector<T> >(3, std::vector<T>(3, T()));

		//.....�v�f�ʐ�.....
		T A = Triangle2DSpace(_nodes[_element[0]], _nodes[_element[1]], _nodes[_element[2]]);

		//.....Ce�}�g���N�X.....
		Ce[0][0] = _rho * _c*A*_t / 6.0;	Ce[0][1] = _rho * _c*A*_t / 12.0;	Ce[0][2] = _rho * _c*A*_t / 12.0;
		Ce[1][0] = _rho * _c*A*_t / 12.0;	Ce[1][1] = _rho * _c*A*_t / 6.0;	Ce[1][2] = _rho * _c*A*_t / 12.0;
		Ce[2][0] = _rho * _c*A*_t / 12.0;	Ce[2][1] = _rho * _c*A*_t / 12.0;	Ce[2][2] = _rho * _c*A*_t / 6.0;
	
		return Ce;
	}
}