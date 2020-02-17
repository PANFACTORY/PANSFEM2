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
	void LinearIsotropicElasticShell(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _v, std::vector<int>& _element, T _E, T _V, T _t) {
		//----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(5*_element.size(), 5*_element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 3);
		for(auto i : _element){
			X = X.Vstack(_x[i].Transpose());
		}
        Vector<T> l = _x[_element[1]] - _x[_element[0]];

		//----------Generate director vector----------
        Matrix<T> v3 = Matrix<T>(0, 3);
        Matrix<T> v1 = Matrix<T>(0, 3);
        Matrix<T> v2 = Matrix<T>(0, 3);
        for(auto i : _element){
            Vector<T> v3i = _v[i].Normal();
            v3 = v3.VStack(v3i.Transpose());
            Vector<T> v1i = VectorProduct(l, v3i).Normal();
            v1 = v1.VStack(v1i.Transpose());
            Vector<T> v2i = VectorProduct(v3i, v1i).Normal();
            v2 = v2.VStack(v2i.Transpose());
        }

        //----------Generate D0 matrix----------
        T k = 1.2;
        Matrix<T> D0 = Matrix<T>(5, 5); 
        D0(0, 0) = 1.0; D0(0, 1) = _V;  D0(0, 2) = T();                 D0(0, 3) = T();                 D0(0, 4) = T();
        D0(1, 0) = _V;  D0(1, 1) = 1.0; D0(1, 2) = T();                 D0(1, 3) = T();                 D0(1, 4) = T();
        D0(2, 0) = T(); D0(2, 1) = T(); D0(2, 2) = 0.5*(1.0 - _V)/k;    D0(2, 3) = T();                 D0(2, 4) = T();
        D0(3, 0) = T(); D0(3, 1) = T(); D0(3, 2) = T();                 D0(3, 3) = 0.5*(1.0 - _V)/k;    D0(3, 4) = T();
        D0(4, 0) = T(); D0(4, 1) = T(); D0(4, 2) = T();                 D0(4, 3) = T();                 D0(4, 4) = 0.5*(1.0 - _V)/k;
        D0 *= _E/(1.0 - _V*_V);

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
            Vector<T> N = SF<T>::N(IC<T>::Points[g].Segment(0, 2));
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g].Segment(0, 2));

			//----------Generate Jacobi matrix and derivative----------
            T zeta = IC<T>::Points[g].Segment(2, 1)(0);
			Matrix<T> J = (dNdr*(X + 0.5*_t*zeta*v3)).VStack(0.5*_t*N.Transpose()*v3);
            T detJ = J.Determinant();
            Matrix<T> invJ = J.Inverse();
            Matrix<T> dNdx = invJ.Block(0, 0, 3, 2)*dNdr;
            Matrix<T> dNzdx = dNdx*zeta + invJ.Block(0, 2, 3, 1)*N.Transpose();

            //----------Generate B matrix----------
            Matrix<T> B = Matrix<T>(6, 5*_element.size());
            for(int i = 0; i < _element.size(); i++){
                B(0, 5*i) = dNdx(0, i); B(0, 5*i + 1) = T();        B(0, 5*i + 2) = T();        B(0, 5*i + 3) = 0.5*_t*dNzdx(0, i)*v1(i, 0);                            B(0, 5*i + 4) = 0.5*_t*dNzdx(0, i)*v2(i, 0);
                B(1, 5*i) = T();        B(1, 5*i + 1) = dNdx(1, i); B(1, 5*i + 2) = T();        B(1, 5*i + 3) = 0.5*_t*dNzdx(1, i)*v1(i, 1);                            B(1, 5*i + 4) = 0.5*_t*dNzdx(1, i)*v2(i, 1);
                B(2, 5*i) = T();        B(2, 5*i + 1) = T();        B(2, 5*i + 2) = dNdx(2, i); B(2, 5*i + 3) = 0.5*_t*dNzdx(2, i)*v1(i, 2);                            B(2, 5*i + 4) = 0.5*_t*dNzdx(2, i)*v2(i, 2);
                B(3, 5*i) = dNdx(1, i); B(3, 5*i + 1) = dNdx(0, i); B(3, 5*i + 2) = T();        B(3, 5*i + 3) = 0.5*_t*(dNzdx(1, i)*v1(i, 0) + dNzdx(0, i)*v1(i, 1));   B(3, 5*i + 4) = 0.5*_t*(dNzdx(1, i)*v2(i, 0) + dNzdx(0, i)*v2(i, 1));
                B(4, 5*i) = dNdx(2, i); B(4, 5*i + 1) = T();        B(4, 5*i + 2) = dNdx(0, i); B(4, 5*i + 3) = 0.5*_t*(dNzdx(0, i)*v1(i, 2) + dNzdx(2, i)*v1(i, 0));   B(4, 5*i + 4) = 0.5*_t*(dNzdx(0, i)*v2(i, 2) + dNzdx(2, i)*v2(i, 0));
                B(5, 5*i) = T();        B(5, 5*i + 1) = dNdx(2, i); B(5, 5*i + 2) = dNdx(1, i); B(5, 5*i + 3) = 0.5*_t*(dNzdx(2, i)*v1(i, 1) + dNzdx(1, i)*v1(i, 2));   B(5, 5*i + 4) = 0.5*_t*(dNzdx(2, i)*v2(i, 1) + dNzdx(1, i)*v2(i, 2));
            }

            //----------Genarate D matrix----------
            Matrix<T> D = Matrix<T>(6, 6);
            

            //----------Update element stiffness matrix----------
            _Ke += B.Transpose()*D*B*detJ*IC::Weights[g][0]*IC::Weights[g][1]*IC::Weights[g][2];
		}
	}
}