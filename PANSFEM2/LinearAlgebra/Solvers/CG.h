//*****************************************************************************
//Title		:CG.h
//Author	:Tanabe Yuta
//Date		:2019/09/09
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <numeric>
#include <chrono>


#include "../Models/CSR.h"


//********************ベクトル和(a+beta*b)*******************
template<class T>
inline std::vector<T> add(std::vector<T> &_a, T _beta, std::vector<T> &_b) {
	std::vector<T> v(_a.size());
	auto ai = _a.begin(), bi = _b.begin();
	for (auto &vi : v) {
		vi = (*ai) + _beta * (*bi);
		++ai;
		++bi;
	}
	return v;
}


//********************ベクトル和(a+beta*b+ganma*c)********************
template<class T>
inline std::vector<T> add(std::vector<T> &_a, T _beta, std::vector<T> &_b, T _ganma, std::vector<T> &_c) {
	std::vector<T> v(_a.size());
	auto ai = _a.begin(), bi = _b.begin(), ci = _c.begin();
	for (auto &vi : v) {
		vi = (*ai) + _beta * (*bi) + _ganma * (*ci);
		++ai;
		++bi;
		++ci;
	}
	return v;
}


//********************ベクトル和差(a+beta(b-ganma*c))********************
template<class T>
inline std::vector<T> addsubstract(std::vector<T> &_a, T _beta, std::vector<T> &_b, T _ganma, std::vector<T> &_c) {
	std::vector<T> v(_a.size());
	auto ai = _a.begin(), bi = _b.begin(), ci = _c.begin();
	for (auto &vi : v) {
		vi = (*ai) + _beta * ((*bi) - _ganma * (*ci));
		++ai;
		++bi;
		++ci;
	}
	return v;
}


//********************ベクトル差(a-b)********************
template<class T>
inline std::vector<T> subtract(std::vector<T> _a, std::vector<T> _b) {
	std::vector<T> v(_b.size());
	auto ai = _a.begin(), bi = _b.begin();
	for (auto &vi : v) {
		vi = (*ai) - (*bi);
		++ai;
		++bi;
	}
	return v;
}


//********************ベクトル差(a-beta*b)*******************
template<class T>
inline std::vector<T> subtract(std::vector<T> &_a, T _beta, std::vector<T> &_b) {
	std::vector<T> v(_a.size());
	auto ai = _a.begin(), bi = _b.begin();
	for (auto &vi : v) {
		vi = (*ai) - _beta * (*bi);
		++ai;
		++bi;
	}
	return v;
}


//********************CG法のメインルーチン********************
template<class T>
std::vector<T> CG(CSR<T> &_A, std::vector<T> &_b, int _itrmax, T _eps) {
	//----------初期化----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> pk = rk;
	T bnorm = std::inner_product(_b.begin(), _b.end(), _b.begin(), T());

	//----------反復計算----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A * pk;
		T alpha = std::inner_product(rk.begin(), rk.end(), rk.begin(), T()) / std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
		std::vector<T> xkp1 = add(xk, alpha, pk);
		std::vector<T> rkp1 = subtract(rk, alpha, Apk);
		T beta = std::inner_product(rkp1.begin(), rkp1.end(), rkp1.begin(), T()) / std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
		std::vector<T> pkp1 = add(rkp1, beta, pk);

		//----------値の更新----------
		xk = xkp1;
		rk = rkp1;
		pk = pkp1;

		//----------収束判定----------
		T rnorm = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************BiCGSTAB法のメインルーチン********************
template<class T>
std::vector<T> BiCGSTAB(CSR<T> &_A, std::vector<T> &_b, int _itrmax, T _eps) {
	//----------初期化----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> rdash = rk;
	std::vector<T> pk = rk;
	T bnorm = std::inner_product(_b.begin(), _b.end(), _b.begin(), T());

	//----------反復計算----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A * pk;
		T rdashdotrk = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
		T alpha = rdashdotrk / std::inner_product(rdash.begin(), rdash.end(), Apk.begin(), T());
		std::vector<T> sk = subtract(rk, alpha, Apk);
		std::vector<T> Ask = _A * sk;
		T omega = std::inner_product(Ask.begin(), Ask.end(), sk.begin(), T()) / std::inner_product(Ask.begin(), Ask.end(), Ask.begin(), T());
		std::vector<T> xkp1 = add(xk, alpha, pk, omega, sk);
		std::vector<T> rkp1 = subtract(sk, omega, Ask);
		T beta = alpha / omega * std::inner_product(rdash.begin(), rdash.end(), rkp1.begin(), T()) / rdashdotrk;
		std::vector<T> pkp1 = addsubstract(rk, beta, pk, omega, Apk);

		//----------xk，rk, pkの更新----------
		xk = xkp1;
		rk = rkp1;
		pk = pkp1;

		//----------収束判定----------
		T rnorm = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************不完全LU(0)分解********************
template<class T>
CSR<T> ILU0(CSR<T> &_A) {
	LILCSR<T> q = LILCSR<T>(_A);

	for (int i = 0, iend = _A.ROWS; i < iend; i++) {
		for (int n = _A.indptr[i], nend = _A.indptr[i + 1]; n < nend; n++) {
			int j = _A.indices[n];
			T qij = _A.data[n];
			for (int k = 0, kend = _A.ROWS; k < kend; k++) {
				if (std::binary_search(_A.indices.begin() + _A.indptr[i], _A.indices.begin() + _A.indptr[i + 1], k)) {
					if (std::binary_search(_A.indices.begin() + _A.indptr[k], _A.indices.begin() + _A.indptr[k + 1], j)) {
						if ((i <= j && k < i) || (i > j && k < j)) {
							qij -= q.get(i, k)*q.get(k, j);
						}
						else {
							break;
						}
					}
				}
			}
			if (i > j) {
				qij /= q.get(j, j);
			}
			q.set(i, j, qij);
		}
	}

	return CSR<T>(q);
}


//********************ILU(0)前処理*******************
template<class T>
std::vector<T> PreILU0(CSR<T> &_A, std::vector<T> &_b) {
	//前進消去（Ly=bを解く）
	std::vector<T> v = std::vector<T>(_b);
	for (int i = 0; i < _b.size(); i++) {
		for (int k = _A.indptr[i]; k < _A.indptr[i + 1]; k++) {
			if (_A.indices[k] < i) {
				v[i] -= _A.data[k] * v[_A.indices[k]];
			}
			else {
				break;
			}
		}
	}

	//後退代入（Ux=yを解く）
	for (int i = _b.size() - 1; i >= 0; i--) {
		for (int k = _A.indptr[i + 1] - 1; k >= _A.indptr[i]; k--) {
			if (_A.indices[k] > i) {
				v[i] -= _A.data[k] * v[_A.indices[k]];
			}
			else {
				break;
			}
		}
		v[i] /= _A.get(i, i);
	}

	return v;
}


//*******************ILU(0)分解前処理付きBiCGSTAB法*******************
template<class T>
std::vector<T> ILU0BiCGSTAB(CSR<T> &_A, CSR<T> &_M, std::vector<T> &_b, int _itrmax, T _eps) {
	//----------初期化----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> rdash = rk;
	std::vector<T> pk = rk;
	T bnorm = std::inner_product(_b.begin(), _b.end(), _b.begin(), T());

	//----------反復計算----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Mpk = PreILU0(_M, pk);		//前処理

		std::vector<T> AMpk = _A * Mpk;
		T rdashdotrk = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
		T alpha = rdashdotrk / std::inner_product(rdash.begin(), rdash.end(), AMpk.begin(), T());
		std::vector<T> sk = subtract(rk, alpha, AMpk);

		std::vector<T> Msk = PreILU0(_M, sk);		//前処理

		std::vector<T> AMsk = _A * Msk;
		T omega = std::inner_product(AMsk.begin(), AMsk.end(), sk.begin(), T()) / std::inner_product(AMsk.begin(), AMsk.end(), AMsk.begin(), T());
		std::vector<T> xkp1 = add(xk, alpha, Mpk, omega, Msk);
		std::vector<T> rkp1 = subtract(sk, omega, AMsk);
		T beta = alpha / omega * std::inner_product(rdash.begin(), rdash.end(), rkp1.begin(), T()) / rdashdotrk;
		std::vector<T> pkp1 = addsubstract(rk, beta, pk, omega, AMpk);

		//----------xk，rk, pkの更新----------
		xk = xkp1;
		rk = rkp1;
		pk = pkp1;

		//----------収束判定----------
		T rnorm = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************SOR法********************
template<class T>
std::vector<T> SOR(CSR<T> &_A, std::vector<T> &_b, T _w, int _itrmax, T _eps) {
	std::vector<T> x = std::vector<T>(_A.ROWS, T());
	for (int itr = 0; itr < _itrmax; itr++) {
		T error = T();
		for (int i = 0; i < _A.ROWS; i++) {
			T Aii = T();
			T tmp = x[i];
			x[i] = _b[i];
			for (int k = _A.indptr[i]; k < _A.indptr[i + 1]; k++) {
				int j = _A.indices[k];
				if (i != j) {
					x[i] -= _A.data[k] * x[j];
				}
				else {
					Aii = _A.data[k];
				}
			}
			x[i] /= Aii;
			x[i] = tmp + _w * (x[i] - tmp);

			error += fabs((tmp - x[i]) / tmp);
		}
		if (error < _eps) {
			std::cout << "\tConvergence:" << itr << std::endl;
			return x;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return x;
}