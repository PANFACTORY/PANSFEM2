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
	//******************************���ʂЂ��݃��f���̍����}�g���N�X�𐶐�******************************
	template<class T>
	void PlaneStrainTri(std::vector<std::vector<T> >& _Ke, std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _E, T _V, T _t) {
		//.....�v�f�����}�g���N�X�̊m��.....
		_Ke = std::vector<std::vector<double> >(6, std::vector<double>(6, T()));
		
		//.....�v�f�ʐ�.....
		double A = 0.5*((_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0])*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) - (_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1])*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]));

		//.....B�}�g���N�X.....
		double B[3][6];
		B[0][0] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[0][1] = 0.0;																B[0][2] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[0][3] = 0.0;																B[0][4] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;	B[0][5] = 0.0;
		B[1][0] = 0.0;																B[1][1] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[1][2] = 0.0;																B[1][3] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[1][4] = 0.0;																B[1][5] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;
		B[2][0] = 0.5*(_nodes[_element[2]].x[0] - _nodes[_element[1]].x[0]) / A;	B[2][1] = 0.5*(_nodes[_element[1]].x[1] - _nodes[_element[2]].x[1]) / A;	B[2][2] = 0.5*(_nodes[_element[0]].x[0] - _nodes[_element[2]].x[0]) / A;	B[2][3] = 0.5*(_nodes[_element[2]].x[1] - _nodes[_element[0]].x[1]) / A;	B[2][4] = 0.5*(_nodes[_element[1]].x[0] - _nodes[_element[0]].x[0]) / A;	B[2][5] = 0.5*(_nodes[_element[0]].x[1] - _nodes[_element[1]].x[1]) / A;

		//.....D�}�g���N�X.....
		double coef = _E / ((1.0 - 2.0*_V)*(1.0 + _V));
		double D[3][3];
		D[0][0] = coef * (1.0 - _V);	D[0][1] = coef * _V;			D[0][2] = 0.0;
		D[1][0] = D[0][1];				D[1][1] = coef * (1.0 - _V);	D[1][2] = 0.0;
		D[2][0] = D[0][2];				D[2][1] = D[1][2];				D[2][2] = coef * 0.5*(1.0 - 2.0*_V);

		//.....DB�}�g���N�X.....
		double DB[3][6];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				DB[i][j] = 0.0;
				for (int k = 0; k < 3; k++) {
					DB[i][j] += D[i][k] * B[k][j];
				}
			}
		}

		//.....Ke(=��BtDBdv)�}�g���N�X.....
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