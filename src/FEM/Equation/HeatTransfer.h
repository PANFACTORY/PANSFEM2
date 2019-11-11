//*****************************************************************************
//Title		:src/FEM/Equation/HeatTransfer.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Heat transfer matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void HeatTransfer(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _alpha, T _t) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(_element.size(), _element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 2; j++){
				X(i, j) = _x[_element[i]](j);
			}
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();
			Matrix<T> B = dXdr.Inverse() * dNdr;

			//----------Update element stiffness matrix----------
			_Ke += B.Transpose()*B*J*_alpha*_t*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}


	//******************************Heat capacity matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void HeatCapacity(Matrix<T>& _Ce, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _rho, T _c, T _t) {
		//----------Get element matrix----------
		_Ce = Matrix<T>(_element.size(), _element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 2; j++){
				X(i, j) = _x[_element[i]](j);
			}
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();

			//----------Make C matrix----------
			_Ce += N*N.Transpose()*J*_rho*_c*_t*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}
}