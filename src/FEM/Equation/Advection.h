//*****************************************************************************
//	Title		:	src/FEM/Equation/Advection.h
//	Author		:	Tanabe Yuta
//	Date		:	2019/10/03
//	Copyright	:	(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Get element advection matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Advection(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _cx, T _cy) {
		assert(_doulist.size() == 1);
		
		_Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> c = Vector<T>({ _cx, _cy });

			_Ke += N*c.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element advection matrix of SUPG******************************
	template<class T, template<class>class SF, template<class>class IC>
	void AdvectionSUPG(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _ax, T _ay, T _k) {
		assert(_doulist.size() == 1);
		
		_Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		Vector<T> a = Vector<T>({ _ax, _ay });
		
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> dNdXa = dNdX.Transpose()*a; 
			T sum = T();
			for(int i = 0; i < _element.size(); i++){
				sum += fabs(dNdXa(i))/a.Norm();
			}
			T he = 2.0/sum;
			T alpha = 0.5*a.Norm()*he/_k;
			T tau = 0.5*he/a.Norm();
			if(alpha <= 3.0){
				tau *= alpha/3.0;
			} else{
				tau *= 1.0;
			}
			
			_Ke += tau*dNdX.Transpose()*a*a.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}		
	}


	//******************************Get element advection matrix of shock capturing******************************
	template<class T, template<class>class SF, template<class>class IC>
	void AdvectionShockCapturing(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _ax, T _ay, T _k) {
		assert(_doulist.size() == 1);
		
		_Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		Vector<T> a = Vector<T>({ _ax, _ay });
		
		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> dNdXa = dNdX.Transpose()*a; 
			T sum = T();
			for(int i = 0; i < _element.size(); i++){
				sum += fabs(dNdXa(i))/a.Norm();
			}
			T he = 2.0/sum;
			T alpha = 0.5*a.Norm()*he/_k;
			T tau = 0.5*a.Norm()*he;
			if(alpha <= 3.0){
				tau *= alpha/3.0;
			} else{
				tau *= 1.0;
			}

			_Ke += tau*dNdX.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}		
	}


	//******************************Get element diffusion matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Diffusion(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _k) {
		assert(_doulist.size() == 1);
		
		_Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			_Ke += _k*dNdX.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element mass matrix******************************
	template<class T, template<class>class SF, template<class>class IC>
	void Mass(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x) {
		assert(_doulist.size() == 1);
		
		_Ce = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();

			_Ce += N*N.Transpose()*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element mass matrix of SUPG******************************
	template<class T, template<class>class SF, template<class>class IC>
	void MassSUPG(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _ax, T _ay, T _k) {
		assert(_doulist.size() == 1);
		
		_Ce = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

		Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);
		}
		
		Vector<T> a = Vector<T>({ _ax, _ay });
		
		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SF<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> dNdXa = dNdX.Transpose()*a; 
			T sum = T();
			for(int i = 0; i < _element.size(); i++){
				sum += fabs(dNdXa(i))/a.Norm();
			}
			T he = 2.0/sum;
			T alpha = 0.5*a.Norm()*he/_k;
			T tau = 0.5*he/a.Norm();
			if(alpha <= 3.0){
				tau *= alpha/3.0;
			} else{
				tau *= 1.0;
			}

			_Ce += tau*dNdX.Transpose()*a*N.Transpose()*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}