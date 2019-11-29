//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/BoundaryCondition.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
	//******************************Set Dirichlet boundary conditions******************************
	template<class T>
	void SetDirichlet(LILCSR<T>& _K, std::vector<T>& _F, std::vector<int>& _isufixed, std::vector<T>& _u, T _alpha) {
		assert(_isufixed.size() == _u.size());

		for (int i = 0; i < _isufixed.size(); i++) {
			T Kii = _K.get(_isufixed[i], _isufixed[i]);
			_F[_isufixed[i]] = _alpha * Kii*_u[i];
			_K.set(_isufixed[i], _isufixed[i], _alpha*Kii);
		}
	}


	//******************************Set Dirichlet boundary conditions******************************
	template<class T>
	void SetDirichlet(LILCSR<T>& _K, LILCSR<T>& _M, std::vector<int>& _isufixed, std::vector<T>& _u, T _alpha) {
		assert(_isufixed.size() == _u.size());

		for (int i = 0; i < _isufixed.size(); i++) {
			_K.set(_isufixed[i], _isufixed[i], _alpha*_K.get(_isufixed[i], _isufixed[i]));
			//_M.set(_isufixed[i], _isufixed[i], _alpha*_M.get(_isufixed[i], _isufixed[i]));
		}
	}


	//******************************Set Neumann boundary conditions******************************
	template<class T>
	void SetNeumann(std::vector<T>& _F, std::vector<int>& _isqfixed, std::vector<T>& _q) {
		assert(_isqfixed.size() == _q.size());

		for (int i = 0; i < _isqfixed.size(); i++) {
			_F[_isqfixed[i]] += _q[i];
		}
	}
}