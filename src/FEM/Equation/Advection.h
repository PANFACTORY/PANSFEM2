//*****************************************************************************
//Title		:src/FEM/Equation/Advection.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Get element advection matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Advection(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _cx, T _cy) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(_element.size(), _element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = (dNdr*X).Transpose();
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Generate advection velocity----------
			Vector<T> c = Vector<T>({ _cx, _cy });

			//----------Update element advection matrix----------
			_Ke += N*c.Transpose()*dNdX*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}


	//******************************Get element advection matrix of SUPG******************************
	template<class T, template<class>class SF, template<class>class IC>
	void AdvectionSUPG(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _ax, T _ay, T _k) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(_element.size(), _element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Generate advection velocity----------
		Vector<T> a = Vector<T>({ _ax, _ay });
		
		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = (dNdr*X).Transpose();
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Get stability parameter----------
			Vector<T> dNdXa = dNdX.Transpose()*(a/a.Norm()); 
			T sum = T();
			for(int i = 0; i < _element.size(); i++){
				sum += fabs(dNdXa(i));
			}
			T he = 2.0/sum;
			T alpha = 0.5*a.Norm()*he/_k;
			T tau = 0.5*he*(1.0/tanh(alpha) - 1.0/alpha)/a.Norm();

			//----------Update element advection matrix----------
			_Ke += tau*dNdX.Transpose()*a*a.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];	
		}		
	}


	//******************************Get element advection matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Diffusion(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _k) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(_element.size(), _element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = (dNdr*X).Transpose();
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Update element advection matrix----------
			_Ke += dNdX.Transpose()*_k*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element mass matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Mass(Matrix<T>& _Ce, std::vector<Vector<T> >& _x, std::vector<int>& _element) {
		//----------Get element matrix----------
		_Ce = Matrix<T>(_element.size(), _element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = (dNdr*X).Transpose();
			T J = dXdr.Determinant();

			//----------Make C matrix----------
			_Ce += N*N.Transpose()*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}