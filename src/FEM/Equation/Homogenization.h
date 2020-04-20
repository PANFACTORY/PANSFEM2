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
	void HomogenizePlaneStrainBodyForce(Matrix<T>& _Fes, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _E, T _V, T _t) {
		assert(_doulist.size() == 2);

		_Fes = Matrix<T>(2*_element.size(), 3);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 2*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 2*i + 1);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		Matrix<T> D = Matrix<T>(3, 3);
		D(0, 0) = 1.0 - _V;	D(0, 1) = _V;		D(0, 2) = T();
		D(1, 0) = D(0, 1);	D(1, 1) = 1.0 - _V;	D(1, 2) = T();
		D(2, 0) = D(0, 2);	D(2, 1) = D(1, 2);	D(2, 2) = 0.5*(1.0 - 2.0*_V);
		D *= _E/((1.0 - 2.0*_V)*(1.0 + _V));

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = dNdX(0, n);	B(0, 2*n + 1) = T();			
				B(1, 2*n) = T();		B(1, 2*n + 1) = dNdX(1, n);	
				B(2, 2*n) = dNdX(1, n);	B(2, 2*n + 1) = dNdX(0, n);	
			}

			_Fes += B.Transpose()*D*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //********************Make homogenized elemental Constitutive********************
    template<class T, template<class>class SF, template<class>class IC>
	Matrix<T> HomogenizePlaneStrainConstitutive(std::vector<Vector<T> >& _x, std::vector<int>& _element, std::vector<Vector<T> >& _chi0, std::vector<Vector<T> >& _chi1, std::vector<Vector<T> >& _chi2, T _E, T _V, T _t) {
		Matrix<T> C = Matrix<T>(3, 3);

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

        Matrix<T> CHI = Matrix<T>(2*_element.size(), 3);
        for(int i = 0; i < _element.size(); i++){
            CHI(2*i, 0) = _chi0[_element[i]](0); 		CHI(2*i, 1) = _chi1[_element[i]](0);		CHI(2*i, 2) = _chi2[_element[i]](0);
            CHI(2*i + 1, 0) = _chi0[_element[i]](1);	CHI(2*i + 1, 1) = _chi1[_element[i]](1);	CHI(2*i + 1, 2) = _chi2[_element[i]](1);
		}

		Matrix<T> D = Matrix<T>(3, 3);
		D(0, 0) = 1.0 - _V;	D(0, 1) = _V;		D(0, 2) = T();
		D(1, 0) = D(0, 1);	D(1, 1) = 1.0 - _V;	D(1, 2) = T();
		D(2, 0) = D(0, 2);	D(2, 1) = D(1, 2);	D(2, 2) = 0.5*(1.0 - 2.0*_V);
		D *= _E/((1.0 - 2.0*_V)*(1.0 + _V));

		Matrix<T> I = Identity<T>(3);

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = dNdX(0, n);	B(0, 2*n + 1) = T();			
				B(1, 2*n) = T();		B(1, 2*n + 1) = dNdX(1, n);	
				B(2, 2*n) = dNdX(1, n);	B(2, 2*n + 1) = dNdX(0, n);	
			}

			C += D*(I - B*CHI)*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}

		return C;
	}


	//********************Check homogenization********************
    template<class T, template<class>class SF, template<class>class IC>
	Matrix<T> HomogenizePlaneStrainCheck(std::vector<Vector<T> >& _x, std::vector<int>& _element, std::vector<Vector<T> >& _chi0, std::vector<Vector<T> >& _chi1, std::vector<Vector<T> >& _chi2, T _t) {
		Matrix<T> C = Matrix<T>(3, 3);

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

        Matrix<T> CHI = Matrix<T>(2*_element.size(), 3);
        for(int i = 0; i < _element.size(); i++){
            CHI(2*i, 0) = _chi0[_element[i]](0); 		CHI(2*i, 1) = _chi1[_element[i]](0);		CHI(2*i, 2) = _chi2[_element[i]](0);
            CHI(2*i + 1, 0) = _chi0[_element[i]](1);	CHI(2*i + 1, 1) = _chi1[_element[i]](1);	CHI(2*i + 1, 2) = _chi2[_element[i]](1);
		}

		Matrix<T> I = Identity<T>(3);

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = dNdX(0, n);	B(0, 2*n + 1) = T();			
				B(1, 2*n) = T();		B(1, 2*n + 1) = dNdX(1, n);	
				B(2, 2*n) = dNdX(1, n);	B(2, 2*n + 1) = dNdX(0, n);	
			}

			C += (I - B*CHI)*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}

		return C;
	}
}