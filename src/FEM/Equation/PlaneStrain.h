//*****************************************************************************
//Title		:src/FEM/Equation/PlaneStrain.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Make element stiffness matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void PlaneStrainStiffness(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _E, T _V, T _t) {
		//----------Initialize element stiffness matrix----------
		_Ke = Matrix<T>(2*_element.size(), 2*_element.size());

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

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
			//----------Get difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse() * dNdr;

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(3, 2 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = dNdX(0, n);	B(0, 2 * n + 1) = 0.0;			
				B(1, 2 * n) = 0.0;			B(1, 2 * n + 1) = dNdX(1, n);	
				B(2, 2 * n) = dNdX(1, n);	B(2, 2 * n + 1) = dNdX(0, n);	
			}

			//----------Update element stiffness matrix----------
			_Ke += B.Transpose()*D*B*J*_t*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}


	//******************************Make element mass matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void PlaneStrainMass(Matrix<T>& _Me, std::vector<Vector<T> >& _x, std::vector<int>& _element, T _rho, T _t) {
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
			Matrix<T> dXdr = dNdr * X;
			T J = dXdr.Determinant();

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(2, 2 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = N(n);	B(0, 2 * n + 1) = 0.0;			
				B(1, 2 * n) = 0.0;	B(1, 2 * n + 1) = N(n);	
			}

			//----------Make C matrix----------
			_Me += B.Transpose()*B*J*_rho*_t*IC<T>::Weights[g][0] * IC<T>::Weights[g][1];
		}
	}


	//******************************Make surface force vector******************************
	template<class T, template<class>class SF, template<class>class IC, class F>
	void PlaneStrainSurfaceForce(Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<int>& _element, F _f, T _t){
		//----------Initialize equivalent nodal force vector----------
		_Fe = Vector<T>(2*_element.size());

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
			Vector<T> x = X.Transpose()*N;
			Matrix<T> dXdr = dNdr*X;
			T dl = sqrt((dXdr*dXdr.Transpose())(0, 0));

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(2, 2 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = N(n);	B(0, 2*n + 1) = 0.0;			
				B(1, 2*n) = 0.0;	B(1, 2*n + 1) = N(n);	
			}

			//----------Get surface force vector----------
			Vector<T> p = _f(x);

			//----------Make Fe Vector----------
			_Fe += B.Transpose()*p*dl*_t*IC<T>::Weights[g][0];
		}
	}


	//******************************Make body force vector******************************
	template<class T, template<class>class SF, template<class>class IC, class F>
	void PlaneStrainBodyForce(Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<int>& _element, F _f, T _t){
		//----------Initialize element mass matrix----------
		_Fe = Vector<T>(2*_element.size());

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
			Vector<T> x = X.Transpose()*N;
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();

			//----------Generate B matrix----------
			Matrix<T> B = Matrix<T>(2, 2 * _element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = N(n);	B(0, 2 * n + 1) = 0.0;			
				B(1, 2 * n) = 0.0;	B(1, 2 * n + 1) = N(n);	
			}

			//----------Get body force vector----------
			Vector<T> p = _f(x);

			//----------Make Fe matrix----------
			_Fe += B.Transpose()*p*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}










	//******************************Make element stiffness matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void PlaneStrainStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _E, T _V, T _t) {
		assert(_doulist.size() == 2);

		_Ke = Matrix<T>(2*_element.size(), 2*_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 2*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 2*i + 1);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

		Matrix<T> D = Matrix<T>(3, 3);
		D(0, 0) = 1.0 - _V;	D(0, 1) = _V;		D(0, 2) = T();
		D(1, 0) = D(0, 1);	D(1, 1) = 1.0 - _V;	D(1, 2) = T();
		D(2, 0) = D(0, 2);	D(2, 1) = D(1, 2);	D(2, 2) = 0.5*(1.0 - 2.0*_V);
		D *= _E / ((1.0 - 2.0*_V)*(1.0 + _V));

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Matrix<T> B = Matrix<T>(3, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = dNdX(0, n);	B(0, 2 * n + 1) = T();			
				B(1, 2 * n) = T();			B(1, 2 * n + 1) = dNdX(1, n);	
				B(2, 2 * n) = dNdX(1, n);	B(2, 2 * n + 1) = dNdX(0, n);	
			}

			_Ke += B.Transpose()*D*B*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Make surface force vector******************************
	template<class T, template<class>class SF, template<class>class IC, class F>
	void PlaneStrainSurfaceForce(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, F _f, T _t){
		assert(_doulist.size() == 2);

		_Fe = Vector<T>(2*_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 2*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 2*i + 1);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Vector<T> x = X.Transpose()*N;
			Matrix<T> dXdr = dNdr*X;
			T dl = sqrt((dXdr*dXdr.Transpose())(0, 0));

			Matrix<T> B = Matrix<T>(2, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2*n) = N(n);	B(0, 2*n + 1) = T();			
				B(1, 2*n) = T();	B(1, 2*n + 1) = N(n);	
			}

			_Fe += B.Transpose()*_f(x)*dl*_t*IC<T>::Weights[g][0];
		}
	}


	//******************************Make body force vector******************************
	template<class T, template<class>class SF, template<class>class IC, class F>
	void PlaneStrainBodyForce(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, F _f, T _t){
		assert(_doulist.size() == 2);

		_Fe = Vector<T>(2*_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 2*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 2*i + 1);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Vector<T> x = X.Transpose()*N;
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();

			Matrix<T> B = Matrix<T>(2, 2*_element.size());
			for (int n = 0; n < _element.size(); n++) {
				B(0, 2 * n) = N(n);	B(0, 2 * n + 1) = T();			
				B(1, 2 * n) = T();	B(1, 2 * n + 1) = N(n);	
			}

			_Fe += B.Transpose()*_f(x)*J*_t*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}