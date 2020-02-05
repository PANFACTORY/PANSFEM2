//*****************************************************************************
//Title		:src/FEM/Equation/Stokes.h
//Author	:Tanabe Yuta
//Date		:2020/02/05
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //******************************Get element advection matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void StokesStatic(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(2*_elementu.size() + _elementp.size(), 2*_elementu.size() + _elementp.size();

		//----------Generate cordinate matrix X of velosity u----------
		Matrix<T> Xu = Matrix<T>(0, 2);
		for(int i = 0; i < _elementu.size(); i++){
			Xu = Xu.Vstack(_x[_elementu[i]].Transpose());
		}

        //----------Generate cordinate matrix X of pressure p----------
		Matrix<T> Xp = Matrix<T>(0, 2);
		for(int i = 0; i < _elementp.size(); i++){
			Xp = Xp.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
            //----------Get shape function for pressure p----------
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dMdr;

			//----------Get shape function for velocity u----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dNdr;

            //----------

			
			//----------Update element advection matrix----------
			_Ke += N*c.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}