//*****************************************************************************
//Title		:LinearAlgebra/LAOperation.h
//Author	:Tanabe Yuta
//Date		:2019/10/07
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <algorithm>
#include <cassert>


namespace PANSFEM2 {
	//***********************************************
	//LinearAlgebra Operation
	//***********************************************


	//********************{a}+{b}********************
	template<class T1, class T2>
	std::vector<T1> operator+(const std::vector<T1>& _a, const std::vector<T2>& _b) {
		assert(_a.size() == _b.size());

		std::vector<T1> c(_a.size());
		std::transform(_a.begin(), _a.end(), _b.begin(), c.begin(), [](T1 _ai, T2 _bi) {return _ai + _bi; });
		return c;
	}


	//********************{a}-{b}********************
	template<class T1, class T2>
	std::vector<T1> operator-(const std::vector<T1>& _a, const std::vector<T2>& _b) {
		assert(_a.size() == _b.size());

		std::vector<T1> c(_a.size());
		std::transform(_a.begin(), _a.end(), _b.begin(), c.begin(), [](T1 _ai, T2 _bi) {return _ai - _bi; });
		return c;
	}


	//********************ƒ¿*{a}********************
	template<class T1, class T2>
	std::vector<T2> operator*(T1 _alpha, const std::vector<T2>& _a) {
		std::vector<T2> c(_a.size());
		std::transform(_a.begin(), _a.end(), c.begin(), [=](T2 _ai) {return _alpha * _ai; });
		return c;
	}


	//********************{a}*ƒ¿********************
	template<class T1, class T2>
	std::vector<T1> operator*(const std::vector<T1>& _a, T2 _alpha) {
		std::vector<T1> c(_a.size());
		std::transform(_a.begin(), _a.end(), c.begin(), [=](T1 _ai) {return _alpha * _ai; });
		return c;
	}


	//********************{a}/ƒ¿********************
	template<class T1, class T2>
	std::vector<T1> operator/(const std::vector<T1>& _a, T2 _alpha) {
		std::vector<T1> c(_a.size());
		std::transform(_a.begin(), _a.end(), c.begin(), [=](T1 _ai) {return _ai / _alpha; });
		return c;
	}


	//********************{a}+={b}********************
	template<class T1, class T2>
	std::vector<T1>& operator+=(std::vector<T1>& _a, const std::vector<T2>& _b) {
		assert(_a.size() == _b.size());

		std::transform(_a.begin(), _a.end(), _b.begin(), _a.begin(), [](T1 _ai, T2 _bi) {return _ai + _bi; });
		return _a;
	}


	//********************{a}-={b}********************
	template<class T1, class T2>
	std::vector<T1>& operator-=(std::vector<T1>& _a, const std::vector<T2>& _b) {
		assert(_a.size() == _b.size());

		std::transform(_a.begin(), _a.end(), _b.begin(), _a.begin(), [](T1 _ai, T2 _bi) {return _ai - _bi; });
		return _a;
	}


	//********************{a}*={b}********************
	template<class T1, class T2>
	std::vector<T1>& operator*=(std::vector<T1>& _a, T2 _alpha) {
		std::transform(_a.begin(), _a.end(), _a.begin(), [=](T1 _ai) {return _ai * _alpha; });
		return _a;
	}


	//********************{a}*={b}********************
	template<class T1, class T2>
	std::vector<T1>& operator/=(std::vector<T1>& _a, T2 _alpha) {
		std::transform(_a.begin(), _a.end(), _a.begin(), [=](T1 _ai) {return _ai / _alpha; });
		return _a;
	}


	//***********************************************
	//Vector Operation
	//***********************************************


	//*******************{a}*{b}********************
	template<class T1, class T2>
	T1 operator*(const std::vector<T1>& _a, const std::vector<T2>& _b) {
		assert(_a.size() == _b.size());
		return std::inner_product(_a.begin(), _a.end(), _b.begin(), T1());
	}


	//*******************<<{a}********************
	template<class T>
	std::ostream& operator<<(std::ostream & _out, const std::vector<T>& _a) {
		for (const auto& ai : _a) {
			_out << ai << std::endl;
		}
		return _out;
	}


	//********************{a}.T********************
	template<class T>
	std::vector<std::vector<T> > Transpose(const std::vector<T>& _a) {
		std::vector<std::vector<T> > c = std::vector<std::vector<T> >(1, std::vector<T>(_a));
		return c;
	}


	//***********************************************
	//Tensor Operation
	//***********************************************


	//********************[A]*{b}********************
	template<class T1, class T2>
	std::vector<T1> operator*(const std::vector<std::vector<T1> >& _A, const std::vector<T2>& _b) {
		std::vector<T1> c(_A.size());
		std::transform(_A.begin(), _A.end(), c.begin(), [&](const std::vector<T1>& _a) {return _a * _b;	});
		return c;
	}


	//********************[A]*[B]********************
	template<class T1, class T2>
	std::vector<std::vector<T1> > operator*(const std::vector<std::vector<T1> >& _A, const std::vector< std::vector<T2> >& _B) {
		std::vector<std::vector<T1> > C = std::vector<std::vector<T1> >(_A.size(), std::vector<T1>(_B[0].size(), T1()));
		for (int i = 0; i < _A.size(); i++) {
			for (int j = 0; j < _B[0].size(); j++) {
				for (int k = 0; k < _A[0].size(); k++) {
					C[i][j] += _A[i][k] * _B[k][j];
				}
			}
		}
		return C;
	}


	//********************<<[A]********************
	template<class T>
	std::ostream& operator<<(std::ostream & _out, const std::vector<std::vector<T> >& _A) {
		for (const auto& Ai : _A) {
			for (const auto& Aij : Ai) {
				_out << Aij << "\t";
			}
			_out << std::endl;
		}
		return _out;
	}


	//********************[A].T********************
	template<class T>
	std::vector<std::vector<T> > Transpose(const std::vector<std::vector<T> >& _A) {
		std::vector<std::vector<T> > C = std::vector<std::vector<T> >(_A[0].size(), std::vector<T>(_A.size()));
		for (int i = 0; i < C.size(); i++) {
			for (int j = 0; j < C[0].size(); j++) {
				C[i][j] = _A[j][i];
			}
		}
		return C;
	}
}