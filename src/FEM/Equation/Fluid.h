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
	template<class T, template<class>class SF, template<class>class IC>
	void FluidMass(Matrix<T>& _Me, std::vector<Vector<T> >& _x, std::vector<int>& _element) {
		//----------Initialize element mass matrix----------
		_Me = Matrix<T>(2*_element.size(), 2*_element.size());

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

			//----------Get Jacobian----------
			T J = (dNdr*X).Determinant();

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(2, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = N(n);	B(0, 2*n + 1) = T();			
				B(1, 2*n) = T();	B(1, 2*n + 1) = N(n);	
			}

			//----------Make Me matrix----------
			_Me += B.Transpose()*B*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //********************Convective term matrix(Gradient form)********************
    template<class T, template<class>class SF, template<class>class IC>
	void Convective(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<int>& _element){
        //----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(2*_element.size(), 2*_element.size());
		
		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

        //----------Generate velocity matrix U----------
		Matrix<T> U = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			U = U.Vstack(_u[_element[i]].Transpose());
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
            Vector<T> u = (U.Transpose())*N;
			Matrix<T> Z = Matrix<T>(2, 4);
			Z(0, 0) = u(0);	Z(0, 1) = T();	Z(0, 2) = u(1);	Z(0, 3) = T();
			Z(1, 0) = T();	Z(1, 1) = u(0);	Z(1, 2) = T();	Z(1, 3) = u(1);

			//----------Generate B and C matrix----------
			Matrix<T> A = Matrix<T>(2, 2*_element.size());
			Matrix<T> B = Matrix<T>(4, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				A(0, 2*n) = N(n);	A(0, 2*n + 1) = T();		
				A(1, 2*n) = T();	A(1, 2*n + 1) = N(n);
				
				B(0, 2*n) = dNdX(0, n);	B(0, 2*n + 1) = T();		
				B(1, 2*n) = T();		B(1, 2*n + 1) = dNdX(0, n);
				B(2, 2*n) = dNdX(1, n);	B(2, 2*n + 1) = T();		
				B(3, 2*n) = T();		B(3, 2*n + 1) = dNdX(1, n);
			}

			//----------Update element stiffness matrix----------
			_Ke += A.Transpose()*Z*B*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1]*IC<T>::Weights[g][2];
		}
    }


    //********************Diffusion term matrix********************
	template<class T, template<class>class SF, template<class>class IC>
	void Diffusion(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(2*_element.size(), 2*_element.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Get B matrix----------
			Matrix<T> B = Matrix<T>(4, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = dNdX(0, n);	B(0, 2*n + 1) = T();		
				B(1, 2*n) = T();		B(1, 2*n + 1) = dNdX(0, n);
				B(2, 2*n) = dNdX(1, n);	B(2, 2*n + 1) = T();		
				B(3, 2*n) = T();		B(3, 2*n + 1) = dNdX(1, n);
			}

			//----------Update element stiffness matrix----------
			_Ke += B.Transpose()*B*J*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}


    //********************Poisson equation for pressure********************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void Poisson(Matrix<T>& _Ke, Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, std::vector<int>& _elementu, std::vector<int>& _elementp) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(_elementp.size(), _elementp.size());
		_Fe = Vector<T>(_elementp.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _elementp.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

		//----------Generate cordinate matrix U----------
		Vector<T> U = Vector<T>();
		for(int i = 0; i < _elementu.size(); i++){
			U = U.Vstack(_u[_element[i]]);
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function for pressure----------
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Update Ke----------
			_Ke += dNdX.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];


			//----------Get shape function for velocity----------
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;
			Matrix<T> B = Matrix<T>(1, 2*_elementu.size());
			for(int n = 0; n < _elementu.size(); n++){
				B(0, 2*n) = dMdX(0, n);	B(0, 2*n + 1) = dMdX(1, n);
			}

			//----------Update Fe----------
			_Fe += N*B*U*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //********************Coefficient matrix for updating velocity********************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void UpdateVelocity(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp) {
		//----------Initialize element matrix----------
		_Ke = Matrix<T>(2*_elementu.size(), _elementp.size());

		//----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _elementp.size(); i++){
			X = X.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function for pressure----------
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);

			//----------Get shape function for velocity----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Get B matrix----------
			Matrix<T> B = Matrix<T>(2, 2*_elementu.size());
			for (int n = 0; n < _elementu.size(); n++) {
				B(0, 2*n) = M(n);	B(0, 2*n + 1) = T();		
				B(1, 2*n) = T();	B(1, 2*n + 1) = M(n);
			}

			//----------Update element stiffness matrix----------
			_Ke += B.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}