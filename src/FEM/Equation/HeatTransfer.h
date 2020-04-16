//*****************************************************************************
//	Title		:	src/FEM/Equation/HeatTransfer.h
//	Author		:	Tanabe Yuta
//	Date		:	2020/04/16
//	Copyright	:	(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Heat transfer matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void HeatTransfer(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _alpha, T _t) {
		assert(_doulist.size() == 1);

		_Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> B = dXdr.Inverse()*dNdr;

			_Ke += B.Transpose()*B*J*_alpha*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Heat capacity matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void HeatCapacity(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _rho, T _c, T _t) {
		assert(_doulist.size() == 1);

		_Ce = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}
		
		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();

			_Ce += N*N.Transpose()*J*_rho*_c*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}