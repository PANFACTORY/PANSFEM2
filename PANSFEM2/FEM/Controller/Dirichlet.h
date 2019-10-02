//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/Dirichlet.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
	//******************************Dirichletã´äEèåèÇÃê›íË******************************
	template<class T>
	void SetDirichlet(LILCSR<T>& _K, std::vector<T>& _F, std::vector<int>& _isufixed, std::vector<T>& _u, T _alpha) {
		assert(_isufixed.size() == _u.size());

		for (int i = 0; i < _isufixed.size(); i++) {
			T Kii = _K.get(_isufixed[i], _isufixed[i]);
			_F[_isufixed[i]] = _alpha * Kii*_u[i];
			_K.set(_isufixed[i], _isufixed[i], _alpha*Kii);
		}
	}
}