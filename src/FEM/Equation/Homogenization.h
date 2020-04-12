//*****************************************************************************
//  Title		:   src/FEM/Equation/Homogenization.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/04/11
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************Make equevalent nodal force for Homogenize********************
	template<class T, template<class>class SF, template<class>class IC>
	void HomogenizePlaneStrainBodyForce(Matrix<T>& _Fes, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _E, T _V, T _t) {
		//----------Initialize element mass matrix----------
		_Fes = Matrix<T>(2*_element.size(), 3);

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Generate D matrix----------
		Matrix<T> D = Matrix<T>(3, 3);
		D(0, 0) = 1.0 - _V;	D(0, 1) = _V;		D(0, 2) = 0.0;
		D(1, 0) = D(0, 1);	D(1, 1) = 1.0 - _V;	D(1, 2) = 0.0;
		D(2, 0) = D(0, 2);	D(2, 1) = D(1, 2);	D(2, 2) = 0.5*(1.0 - 2.0*_V);
		D *= (_E / ((1.0 - 2.0*_V)*(1.0 + _V)));

		//----------Generate "i"th unit vector----------
		Matrix<T> I = Identity<T>(3);
		
		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = dNdX(0, n);	B(0, 2 * n + 1) = 0.0;			
				B(1, 2 * n) = 0.0;			B(1, 2 * n + 1) = dNdX(1, n);	
				B(2, 2 * n) = dNdX(1, n);	B(2, 2 * n + 1) = dNdX(0, n);	
			}

			//----------Make Fe matrix----------
			_Fes += B.Transpose()*D*I*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //********************Make homogenized elemental Constitutive********************
    template<class T, template<class>class SF, template<class>class IC>
	void HomogenizePlaneStrainConstitutive(Matrix<T>& _C, std::vector<Vector<T> >& _x, std::vector<int>& _element, std::vector<Vector<T> >& _chi0, std::vector<Vector<T> >& _chi1, std::vector<Vector<T> >& _chi2, T _E, T _V, T _t) {
		//----------Initialize element mass matrix----------
		_C = Matrix<T>(3, 3);

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

        //----------Generate Characteristic displacement matrix----------
        Matrix<T> CHI = Matrix<T>(0, 3);
        for(int i = 0; i < _element.size(); i++){
            Matrix<T> CHIi = Matrix<T>(2, 3);
            CHIi(0, 0) = _chi0[_element[i]](0); CHIi(0, 1) = _chi1[_element[i]](0); CHIi(0, 2) = _chi2[_element[i]](0);
            CHIi(1, 0) = _chi0[_element[i]](1); CHIi(1, 1) = _chi1[_element[i]](1); CHIi(1, 2) = _chi2[_element[i]](1);
            CHI = CHI.Vstack(CHIi);
		}

		//----------Generate D matrix----------
		Matrix<T> D = Matrix<T>(3, 3);
		D(0, 0) = 1.0 - _V;	D(0, 1) = _V;		D(0, 2) = 0.0;
		D(1, 0) = D(0, 1);	D(1, 1) = 1.0 - _V;	D(1, 2) = 0.0;
		D(2, 0) = D(0, 2);	D(2, 1) = D(1, 2);	D(2, 2) = 0.5*(1.0 - 2.0*_V);
		D *= (_E / ((1.0 - 2.0*_V)*(1.0 + _V)));

        //----------Generate "i"th unit vector----------
		Matrix<T> I = Identity<T>(3);

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = dNdX(0, n);	B(0, 2 * n + 1) = 0.0;			
				B(1, 2 * n) = 0.0;			B(1, 2 * n + 1) = dNdX(1, n);	
				B(2, 2 * n) = dNdX(1, n);	B(2, 2 * n + 1) = dNdX(0, n);	
			}

			//----------Make Fe matrix----------
			_C += D*(I - B*CHI)*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}