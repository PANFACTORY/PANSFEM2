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
			if(fabs(Kii) < 1.0e-10){
				_F[_isufixed[i]] = _alpha*_u[i];
				_K.set(_isufixed[i], _isufixed[i], _alpha);
			} else{
				_F[_isufixed[i]] = _alpha*Kii*_u[i];
				_K.set(_isufixed[i], _isufixed[i], _alpha*Kii);
			}
		}
	}


	//******************************Set Dirichlet boundary conditions******************************
	template<class T>
	void SetDirichlet(LILCSR<T>& _K, std::vector<int>& _isufixed, std::vector<T>& _u, T _alpha) {
		assert(_isufixed.size() == _u.size());

		for (int i = 0; i < _isufixed.size(); i++) {
			_K.set(_isufixed[i], _isufixed[i], _alpha*_K.get(_isufixed[i], _isufixed[i]));
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


	//******************************Set Periodic boundary conditions******************************
	template<class T>
	void SetPeriodic(LILCSR<T>& _K, std::vector<int>& _ismasterfixed, std::vector<int>& _isslavefixed, T _alpha) {
		assert(_ismasterfixed.size() == _isslavefixed.size());

		for(int i = 0; i < _ismasterfixed.size(); i++) {
			_K.set(_ismasterfixed[i], _ismasterfixed[i], _K.get(_ismasterfixed[i], _ismasterfixed[i]) + _alpha);
			_K.set(_ismasterfixed[i], _isslavefixed[i], -_alpha);
			_K.set(_isslavefixed[i], _ismasterfixed[i], -_alpha);
			_K.set(_isslavefixed[i], _isslavefixed[i], _K.get(_isslavefixed[i], _isslavefixed[i]) + _alpha);
		}
	}


	//******************************Set Initial boundary conditions*******************************
	template<class T>
	void SetInitial(std::vector<Vector<T> >& _values, const std::vector<int>& _field, const std::vector<int>& _isufixed, const std::vector<T>& _ufixed) {
		assert(_isufixed.size() == _ufixed.size() && _values.size() == _field.size() - 1);

		std::vector<T> result = std::vector<T>(_field[_field.size() - 1], T());
		for(int i = 0; i < _isufixed.size(); i++) {
			result[_isufixed[i]] += _ufixed[i];
		}

		int resultindex = 0;
		for (int i = 0; i < _field.size() - 1; i++) {
			std::vector<T> value;
			for (int k = _field[i]; k < _field[i + 1]; k++) {
				value.push_back(result[resultindex]);
				resultindex++;
			}
			_values[i] = Vector<T>(value);
		}
	}
}