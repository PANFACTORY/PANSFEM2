//*****************************************************************************
//  Title       :src/FEM/Equation/Fluid.h
//  Author      :Tanabe Yuta
//  Date        :2019/12/20
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************




#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Mass term matrix(Lumped mass)********************

    //********************Convective term matrix(Gradient form)********************
    template<class T, template<class>class SF, template<class>class IC>
	void Convective_GradientForm(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<int>& _element){
        //----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(2*_element.size(), 2*_element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 2; j++){
				X(i, j) = _x[_element[i]](j);
			}
		}

        //----------Generate velocity matrix U----------
		Matrix<T> U = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			for(int j = 0; j < 2; j++){
				U(i, j) = _u[_element[i]](j);
			}
		}

        //----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

            //----------Get velocity gradient----------
            Matrix<T> Z = (dNdX*U).Transpose();

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(2, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = N(n);	B(0, 2*n + 1) = T();		
				B(1, 2*n) = T();	B(1, 2*n + 1) = N(n);
			}

			//----------Update element stiffness matrix----------
			_Ke += B.Transpose()*Z*B*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1] * IC<T>::Weights[g][2];
		}
    }

    //********************Diffusion term matrix********************

    //********************Coefficient matrix of Poisson equation for pressure********************

    //********************Coefficient vector of Poisson equation for pressure********************

    //********************Coefficient matrix for updating velocity********************
}