//*****************************************************************************
//Title		:PANSFEM2/LinearAlgebra/Models/CSR.h
//Author	:Tanabe Yuta
//Date		:2019/10/01
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <omp.h>


#include "LILCSR.h"


template<class T>
class LILCSR;


template<class T>
class CSR
{
public:
	CSR();
	~CSR();
	CSR(int _rows, int _cols);	//_rows�F�s���C_cols�F��
	CSR(LILCSR<T> _matrix);		//LILCSR����CSR�𐶐�


	const int ROWS;				//�s��
	const int COLS;				//��


	const std::vector<T> operator*(const std::vector<T> &_vec);					//�x�N�g���Ƃ̐�


	template<class T1, class T2>
	friend const CSR<T1> operator+(const CSR<T1>& _m1, const CSR<T2>& _m2);		//�s��Ƃ̘a
	template<class T1, class T2>
	friend const CSR<T1> operator-(const CSR<T1>& _m1, const CSR<T2>& _m2);		//�s��Ƃ̍�
	template<class T1, class T2>
	friend const CSR<T2> operator*(T1 _a, const CSR<T2>& _m);					//�����Ƃ̐�
	template<class T1, class T2>
	friend const CSR<T1> operator*(const CSR<T1>& _m, T2 _a);					//�����Ƃ̐�
	template<class T1, class T2>
	friend const CSR<T1> operator/(const CSR<T1>& _m, T2 _a);					//�����Ƃ̏�
	template<class F>
	friend std::ostream& operator << (std::ostream &_out, const CSR<F> &_mat);	//stream�ɏo��


	bool set(int _row, int _col, T _data);		//�l�̃Z�b�g
	T get(int _row, int _col) const;			//�l�̎擾


	template<class F>
	friend CSR<F> ILU0(CSR<F>& _A);				//�s���SLU(0)����
	template<class F>
	friend std::vector<F> PreILU0(CSR<F> &_A, std::vector<F> &_b);		//�s���SLU(0)����O����
	template<class F>
	friend std::vector<F> SOR(CSR<F> &_A, std::vector<F> &_b, F _w, int _itrmax, F _eps);	//SOR�@


	template<class F>
	friend class LILCSR;


private:
	std::vector<int> indptr;	//indptr[i]<=j<indptr[i+1]��data[j],indices[j]��i�s�ڂ̗v�f�ł���i�T�C�Y�͍s��+1�j
	std::vector<int> indices;	//indices[j]��data[j]������ڂ̗v�f���������idata�̃T�C�Y�ɓ������j
	std::vector<T> data;		//�i�[�f�[�^�i��[���v�f�̐��ɓ������j
};


template<class T>
inline CSR<T>::CSR() : ROWS(0), COLS(0) {}


template<class T>
inline CSR<T>::~CSR() {}


template<class T>
inline CSR<T>::CSR(int _rows, int _cols) : ROWS(_rows), COLS(_cols) {
	this->indptr = std::vector<int>(this->ROWS + 1, 0);
}


template<class T>
inline CSR<T>::CSR(LILCSR<T> _matrix) : ROWS(_matrix.ROWS), COLS(_matrix.COLS) {
	this->indptr = std::vector<int>(this->ROWS + 1, 0);

	for (int i = 0; i < this->ROWS; i++) {
		this->indptr[i + 1] = this->indptr[i] + _matrix.data[i].size();

		std::sort(_matrix.data[i].begin(), _matrix.data[i].end());
		for (auto matrixj : _matrix.data[i]) {
			this->indices.push_back(matrixj.first);
			this->data.push_back(matrixj.second);
		}
	}
}


template<class T>
inline const std::vector<T> CSR<T>::operator*(const std::vector<T> &_vec) {
	std::vector<T> v(this->ROWS, T());

	int iend = this->ROWS;

#pragma omp parallel for
	for (int i = 0; i < iend; ++i) {
		for (int j = this->indptr[i], jend = this->indptr[i + 1]; j < jend; ++j) {
			v[i] += this->data[j] * _vec[this->indices[j]];
		}
	}

	return v;
}


template<class T>
inline bool CSR<T>::set(int _row, int _col, T _data) {
	auto colbegin = this->indices.begin() + this->indptr[_row], colend = this->indices.begin() + this->indptr[_row + 1];
	auto colnow = std::lower_bound(colbegin, colend, _col);
	if (colnow == colend) {
		this->data.insert(this->data.begin() + std::distance(this->indices.begin(), colnow), _data);
		this->indices.insert(colnow, _col);
		for (auto i = this->indptr.begin() + (_row + 1), iend = this->indptr.end(); i != iend; ++i) {
			*i += 1;
		}
		return false;
	}
	else if (*colnow == _col) {
		*(this->data.begin() + std::distance(this->indices.begin(), colnow)) = _data;
		return true;
	}
	else {
		this->data.insert(this->data.begin() + std::distance(this->indices.begin(), colnow), _data);
		this->indices.insert(colnow, _col);
		for (auto i = this->indptr.begin() + (_row + 1), iend = this->indptr.end(); i != iend; ++i) {
			*i += 1;
		}
		return false;
	}

	return false;
}


template<class T>
inline T CSR<T>::get(int _row, int _col) const {
	auto colbegin = this->indices.begin() + this->indptr[_row], colend = this->indices.begin() + this->indptr[_row + 1];
	auto colnow = std::lower_bound(colbegin, colend, _col);
	if (colnow == colend) {
		return T();
	}
	else if (*colnow == _col) {
		return *(this->data.begin() + std::distance(this->indices.begin(), colnow));
	}
	else {
		return T();
	}
}


template<class F>
inline CSR<F> operator*(F _a, const CSR<F>& _m) {
	CSR<F> m = CSR<F>(_m);
	for (auto& datai : m.data) {
		datai *= _a;
	}
	return m;
}


template<class T1, class T2>
inline const CSR<T1> operator+(const CSR<T1>& _m1, const CSR<T2>& _m2) {
	assert(_m1.ROWS == _m2.ROWS && _m1.COLS == _m2.COLS);

	CSR<T1> m = CSR<T1>(_m1);

	for (int i = 0; i < _m2.ROWS; i++) {
		for (int k = _m2.indptr[i]; k < _m2.indptr[i + 1]; k++) {
			int j = _m2.indices[k];
			m.set(i, j, m.get(i, j) + _m2.data[k]);
		}
	}

	return m;
}


template<class T1, class T2>
inline const CSR<T1> operator-(const CSR<T1>& _m1, const CSR<T2>& _m2) {
	assert(_m1.ROWS == _m2.ROWS && _m1.COLS == _m2.COLS);

	CSR<T1> m = CSR<T1>(_m1);

	for (int i = 0; i < _m2.ROWS; i++) {
		for (int k = _m2.indptr[i]; k < _m2.indptr[i + 1]; k++) {
			int j = _m2.indices[k];
			m.set(i, j, m.get(i, j) - _m2.data[k]);
		}
	}

	return m;
}


template<class T1, class T2>
inline const CSR<T2> operator*(T1 _a, const CSR<T2>& _m) {
	CSR<T2> m = CSR<T2>(_m);
	for (auto& datai : m.data) {
		datai *= _a;
	}
	return m;
}


template<class T1, class T2>
inline const CSR<T1> operator*(const CSR<T1>& _m, T2 _a) {
	CSR<T1> m = CSR<T1>(_m);
	for (auto& datai : m.data) {
		datai *= _a;
	}
	return m;
}


template<class T1, class T2>
inline const CSR<T1> operator/(const CSR<T1>& _m, T2 _a) {
	CSR<T1> m = CSR<T1>(_m);
	for (auto& datai : m.data) {
		datai /= _a;
	}
	return m;
}


template<class F>
inline std::ostream & operator<<(std::ostream & _out, const CSR<F>& _mat) {
	for (int i = 0; i < _mat.ROWS; i++) {
		for (int j = 0; j < _mat.COLS; j++) {
			_out << _mat.get(i, j) << "\t";
		}
		_out << std::endl;
	}

	return _out;
}