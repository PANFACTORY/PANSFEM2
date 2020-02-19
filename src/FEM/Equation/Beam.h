//*****************************************************************************
//Title		:src/FEM/Equation/Beam.h
//Author	:Tanabe Yuta
//Date		:2020/02/19
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Linear Isotropic Elastic Beam********************
	template<class T, template<class>class SF, template<class>class IC0, template<class>class IC12>
	void LinearIsotropicElasticBeam(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, std::vector<Vector<T> >& _v1, std::vector<Vector<T> >& _v2, T _E, T _V, T _a, T _b) {
        //----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(6*_element.size(), 6*_element.size());
		
       //----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 3);
		for(auto i : _element){
			X = X.Vstack(_x[i].Transpose());
		}

        //----------Generate director vector----------
        Matrix<T> v1 = Matrix<T>(0, 3);
        Matrix<T> v2 = Matrix<T>(0, 3);
        for(auto i : _element){
            v1 = v1.Vstack(_v1[i].Normal().Transpose());
            v2 = v2.Vstack(_v2[i].Normal().Transpose());
        }

        //----------Generate local D matrix----------
        T k = 1.2;
        Matrix<T> D0 = Matrix<T>(3, 3);
        D0(0, 0) = _E;  D0(0, 1) = T();                 D0(0, 2) = T();
        D0(1, 0) = T(); D0(1, 1) = 0.5*_E/(1.0 + _V)/k; D0(1, 2) = T();
        D0(2, 0) = T(); D0(2, 1) = T();                 D0(2, 2) = 0.5*_E/(1.0 + _V)/k;

        //----------Loop of Gauss Integration----------
        for(int g = 0; g < IC0<T>::N; g++){
            for(int h = 0; h < IC12<T>::N; h++){
                //----------Get difference of shape function----------
                Vector<T> N = SF<T>::N(IC0<T>::Points[g]);
                Matrix<T> dNdr = SF<T>::dNdr(IC0<T>::Points[g]);
                
                //----------Generate Jacobi matrix and derivative----------
                
                //----------Generate B matrix----------
                Matrix<T> B = Matrix<T>(6, 6*_element.size());

                //----------Generate global D matrix----------

                //----------Update element stiffness matrix----------
            }
        } 
    }
}