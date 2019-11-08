//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/Assembling.h
//Author	:Tanabe Yuta
//Date		:2019/10/04
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/LILCSR.h"
#include "../../LinearAlgebra/Models/Matrix.h"


namespace PANSFEM2 {
	//********************Assembling K Matrix from Ke Matrix********************
	template<class T>
	void Assembling(LILCSR<T>& _K, const std::vector<std::vector<T> >& _Ke, const std::vector<int>& _element, const std::vector<int>& _field) {
		int ei = 0;
		for (auto ni : _element) {
			for (int si = _field[ni]; si < _field[ni + 1]; si++) {
				int ej = 0;
				for (auto nj : _element) {
					for (int sj = _field[nj]; sj < _field[nj + 1]; sj++) {
						_K.set(si, sj, _K.get(si, sj) + _Ke[ei][ej]);
						ej++;
					}
				}
				ei++;
			}
		}
	}


	template<class T>
	void Assembling(LILCSR<T>& _K, Matrix<T>& _Ke, const std::vector<int>& _element, const std::vector<int>& _field) {
		int ei = 0;
		for (auto ni : _element) {
			for (int si = _field[ni]; si < _field[ni + 1]; si++) {
				int ej = 0;
				for (auto nj : _element) {
					for (int sj = _field[nj]; sj < _field[nj + 1]; sj++) {
						_K.set(si, sj, _K.get(si, sj) + _Ke(ei, ej));
						ej++;
					}
				}
				ei++;
			}
		}
	}


	//********************Assembling K & F from Ke & Fe********************
	template<class T>
	void Assembling(LILCSR<T>& _K, std::vector<T>& _F, const std::vector<std::vector<T> >& _Ke, const std::vector<T>& _Fe, const std::vector<int>& _element, const std::vector<int>& _field) {
		int ei = 0;
		for (auto ni : _element) {
			for (int si = _field[ni]; si < _field[ni + 1]; si++) {
				int ej = 0;
				for (auto nj : _element) {
					for (int sj = _field[nj]; sj < _field[nj + 1]; sj++) {
						_K.set(si, sj, _K.get(si, sj) + _Ke[ei][ej]);
						ej++;
					}
				}
				_F[si] += _Fe[ei];
				ei++;
			}
		}
	}


	//********************FieldResultToNodeVector********************
	template<class T>
	void FieldResultToNodeValue(const std::vector<T>& _result, std::vector<std::vector<T> >& _values, const std::vector<int>& _field) {
		int resultindex = 0;
		for (int i = 0; i < _field.size() - 1; i++) {
			std::vector<T> value;
			for (int k = _field[i]; k < _field[i + 1]; k++) {
				value.push_back(_result[resultindex]);
				resultindex++;
			}
			_values.push_back(value);
		}
	}


	//********************FieldResultToNodeScaler********************
	template<class T>
	void FieldResultToNodeValue(const std::vector<T>& _result, std::vector<T>& _values, const std::vector<int>& _field) {
		int resultindex = 0;
		for (int i = 0; i < _field.size() - 1; i++) {
			T value = T();
			for (int k = _field[i]; k < _field[i + 1]; k++) {
				value = _result[resultindex];
				resultindex++;
			}
			_values.push_back(value);
		}
	}
}