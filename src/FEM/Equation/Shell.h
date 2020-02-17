//*****************************************************************************
//Title		:src/FEM/Equation/Shell.h
//Author	:Tanabe Yuta
//Date		:2020/02/17
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Linear Isotropic Elastic Shell********************
	template<class T, template<class>class SF, template<class>class IC>
	void LinearIsotropicElasticShell(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _E, T _V, T _t) {
		//----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(5*_element.size(), 5*_element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 3);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}
        Vector<T> l = _x[_element[1]] - _x[_element[0]];

		//----------Generate director vector----------
        Matrix<T> v3 = Matrix<T>(0, 3);
        Matrix<T> v1 = Matrix<T>(0, 3);
        Matrix<T> v2 = Matrix<T>(0, 3);
        for(int i = 0; i < _element.size(); i++){
            Vector<T> v3i;
            v3 = v3.VStack(v3i.Transpose());
            Vector<T> v1i = VectorProduct(l, v3i).Normal();
            v1 = v1.VStack(v1i.Transpose());
            Vector<T> v2i = VectorProduct(v3i, v1i).Normal();
            v2 = v2.VStack(v2i.Transpose());
        }

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
            Vector<T> N = SF<T>::N(IC<T>::Points[g].Segment(0, 2));
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g].Segment(0, 2));

			//----------Generate Jacobi matrix and derivative----------
			Matrix<T> J = (dNdr*(X + 0.5*_t*IC<T>::Points[g].Segment(2, 1)*v3)).VStack(0.5*_t*N.Transpose()*v3);
            T detJ = J.Determinant();

            //----------Generate B matrix----------
            Matrix<T> B = Matrix<T>(5, 5*_element.size());
            for(int i = 0; i < _element.size(); i++){

            }

            //----------Genarate D matrix----------
            Matrix<T> D = Matrix<T>(5, 5);

            //----------Update element stiffness matrix----------
            _Ke += B.Transpose()*D*B*detJ*IC::Weights[g][0]*IC::Weights[g][1]*IC::Weights[g][2];
		}
	}
}