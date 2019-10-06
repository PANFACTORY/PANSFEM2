#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>


#include "CSR.h"


template<class T>
class CSR;


template<class T>
class LILCSR
{
public:
	LILCSR();
	~LILCSR();
	LILCSR(int _rows, int _cols);
	LILCSR(CSR<T> _matrix);						//CSRからLILCSRを生成


	const int ROWS;			//行数
	const int COLS;			//列数


	template<class T1, class T2>
	friend const std::vector<T1> operator*(const LILCSR<T1>& _m, const std::vector<T2> &_vec);		//ベクトルとの積
	template<class T1, class T2>
	friend const LILCSR<T1> operator+(const LILCSR<T1>& _m1, const LILCSR<T2>& _m2);	//行列との和
	template<class T1, class T2>
	friend const LILCSR<T1> operator-(const LILCSR<T1>& _m1, const LILCSR<T2>& _m2);	//行列との差
	template<class T1, class T2>
	friend const LILCSR<T2> operator*(T1 _a, const LILCSR<T2>& _m);	//実数との積
	template<class T1, class T2>
	friend const LILCSR<T1> operator*(const LILCSR<T1>& _m, T2 _a);	//実数との積
	template<class T1, class T2>
	friend const LILCSR<T1> operator/(const LILCSR<T1>& _m, T2 _a);	//実数との商
	template<class F>
	friend std::ostream& operator << (std::ostream &_out, const LILCSR<F> &_mat);		//streamに出力


	bool set(int _row, int _col, T _data);		//値のセット
	T get(int _row, int _col) const;			//値の取得


	template<class F>
	friend class CSR;


private:
	std::vector<std::vector<std::pair<int, T> > > data;
};


template<class T>
inline LILCSR<T>::LILCSR() : ROWS(0), COLS(0) {}


template<class T>
inline LILCSR<T>::~LILCSR() {}


template<class T>
inline LILCSR<T>::LILCSR(int _rows, int _cols) : ROWS(_rows), COLS(_cols) {
	this->data = std::vector<std::vector<std::pair<int, T> > >(_rows);
}


template<class T>
inline LILCSR<T>::LILCSR(CSR<T> _matrix) : ROWS(_matrix.ROWS), COLS(_matrix.COLS) {
	this->data = std::vector<std::vector<std::pair<int, T> > >(_matrix.ROWS);
	for (int i = 0; i < _matrix.ROWS; i++) {
		for (int k = _matrix.indptr[i]; k < _matrix.indptr[i + 1]; k++) {
			this->data[i].push_back(std::pair<int, T>(_matrix.indices[k], _matrix.data[k]));
		}
	}
}


template<class T>
inline bool LILCSR<T>::set(int _row, int _col, T _data) {
	for (auto& dataj : this->data[_row]) {
		if (dataj.first == _col) {
			dataj.second = _data;
			return true;
		}
	}

	this->data[_row].push_back(std::pair<int, T>(_col, _data));
	return false;
}


template<class T>
inline T LILCSR<T>::get(int _row, int _col) const {
	for (const auto& dataj : this->data[_row]) {
		if (dataj.first == _col) {
			return dataj.second;
		}
	}

	return T();
}


template<class T1, class T2>
inline const std::vector<T1> operator*(const LILCSR<T1>& _m, const std::vector<T2>& _vec) {
	assert(_m.COLS == _vec.size());

	std::vector<T1> v(_m.ROWS, T1());

	//#pragma omp parallel for 
	for (int i = 0; i < _m.ROWS; i++) {
		for (auto dataj : _m.data[i]) {
			v[i] += dataj.second * _vec[dataj.first];
		}
	}

	return v;
}


template<class T1, class T2>
inline const LILCSR<T1> operator+(const LILCSR<T1>& _m1, const LILCSR<T2>& _m2) {
	assert(_m1.ROWS == _m2.ROWS && _m1.COLS == _m2.COLS);

	LILCSR<T1> m = LILCSR<T1>(_m1);
	for (int i = 0; i < _m2.ROWS; i++) {
		for (auto dataij : _m2.data[i]) {
			int j = dataij.first;
			m.set(i, j, m.get(i, j) + dataij.second);
		}
	}

	return m;
}


template<class T1, class T2>
inline const LILCSR<T1> operator-(const LILCSR<T1>& _m1, const LILCSR<T2>& _m2) {
	assert(_m1.ROWS == _m2.ROWS && _m1.COLS == _m2.COLS);

	LILCSR<T1> m = LILCSR<T1>(_m1);
	for (int i = 0; i < _m2.ROWS; i++) {
		for (auto dataij : _m2.data[i]) {
			int j = dataij.first;
			m.set(i, j, m.get(i, j) - dataij.second);
		}
	}

	return m;
}


template<class T1, class T2>
inline const LILCSR<T2> operator*(T1 _a, const LILCSR<T2>& _m) {
	LILCSR<T2> m = LILCSR<T2>(_m);
	for (auto& row : m.data) {
		for (auto& col : row) {
			col.second *= _a;
		}
	}
	return m;
}


template<class T1, class T2>
inline const LILCSR<T1> operator*(const LILCSR<T1>& _m, T2 _a) {
	LILCSR<T1> m = LILCSR<T1>(_m);
	for (auto& row : m.data) {
		for (auto& col : row) {
			col.second *= _a;
		}
	}
	return m;
}


template<class T1, class T2>
inline const LILCSR<T1> operator/(const LILCSR<T1>& _m, T2 _a) {
	LILCSR<T1> m = LILCSR<T1>(_m);
	for (auto& row : m.data) {
		for (auto& col : row) {
			col.second /= _a;
		}
	}
	return m;
}


template<class F>
inline std::ostream & operator<<(std::ostream & _out, const LILCSR<F>& _mat) {
	for (int i = 0; i < _mat.ROWS; i++) {
		for (int j = 0; j < _mat.COLS; j++) {
			_out << _mat.get(i, j) << "\t";
		}
		_out << std::endl;
	}

	return _out;
}