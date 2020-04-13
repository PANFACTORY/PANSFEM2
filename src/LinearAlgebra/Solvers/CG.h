//*****************************************************************************
//Title		:LinearAlgebra/Solvers/CG.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <numeric>
#include <chrono>
#include "../Models/CSR.h"


//********************{a}-{b}********************
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


//********************{x}={x}+a{y}********************
template<class T>
inline void xexpay(std::vector<T>& _x, T _a, const std::vector<T>& _y) {
	auto yi = _y.begin();
	for(auto& xi : _x) {
		xi = xi + _a*(*yi);
		++yi;
	}
}


//********************{x}=a{x}+{y}********************
template<class T>
inline void xeaxpy(T _a, std::vector<T>& _x, const std::vector<T>& _y) {
	auto yi = _y.begin();
	for(auto& xi : _x) {
		xi = _a*xi + (*yi);
		++yi;
	}
}


//********************{z}=a{w}+b({x}-{y}+c{z})********************
template<class T>
inline void zeawpbxmypcz(T _a, const std::vector<T>& _w, T _b, const std::vector<T>& _x, const std::vector<T>& _y, T _c, std::vector<T>& _z) {
	auto wi = _w.begin(), xi = _x.begin(), yi = _y.begin();
	for(auto& zi : _z) {
		zi = _a*(*wi) + _b*((*xi) - (*yi) + _c*zi);
		++wi;
		++xi;
		++yi;
	}
}


//********************{x}=a{x}+b{y}+c{z}********************
template<class T>
inline void xeaxpbypcz(T _a, std::vector<T>& _x, T _b, const std::vector<T>& _y, T _c, const std::vector<T>& _z) {
	auto yi = _y.begin(), zi = _z.begin();
	for(auto& xi : _x) {
		xi = _a*xi + _b*(*yi) + _c*(*zi);
		++yi;
		++zi;
	}
}


//********************{z}=a{x}+b{y}********************
template<class T>
inline std::vector<T> zeaxpby(T _a, const std::vector<T>& _x, T _b, const std::vector<T>& _y) {
	std::vector<T> z = std::vector<T>(_x.size());
	auto xi = _x.begin(), yi = _y.begin();
	for(auto& zi : z) {
		zi = _a*(*xi) + _b*(*yi);
		++xi;
		++yi;
	} 
	return z;
}


//********************{z}=a{w}+b{x}+c{y}********************
template<class T>
inline std::vector<T> zeawpbxpcy(T _a, const std::vector<T>& _w, T _b, const std::vector<T>& _x, T _c, const std::vector<T>& _y) {
	std::vector<T> z = std::vector<T>(_x.size());
	auto wi = _w.begin(), xi = _x.begin(), yi = _y.begin();
	for(auto& zi : z) {
		zi = _a*(*wi) + _b*(*xi) + _c*(*yi);
		++wi;
		++xi;
		++yi;
	} 
	return z;
}


//********************{z}=a{v}+b{w}+c{x}+d{y}********************
template<class T>
inline std::vector<T> zeavpbwpcxpdy(T _a, const std::vector<T>& _v, T _b, const std::vector<T>& _w, T _c, const std::vector<T>& _x, T _d, const std::vector<T>& _y) {
	std::vector<T> z = std::vector<T>(_x.size());
	auto vi = _v.begin(), wi = _w.begin(), xi = _x.begin(), yi = _y.begin();
	for(auto& zi : z) {
		zi = _a*(*vi) + _b*(*wi) + _c*(*xi) + _d*(*yi);
		++vi;
		++wi;
		++xi;
		++yi;
	} 
	return z;
}


//********************CG method********************
template<class T>
std::vector<T> CG(CSR<T>& _A, const std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Initialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> pk = rk;
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));
	T rkrk = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());

	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A*pk;
		T alpha = rkrk/std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
		xexpay(xk, alpha, pk);
		xexpay(rk, -alpha, Apk);
		T rkp1rkp1 = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
		T beta = rkp1rkp1/rkrk;
		xeaxpy(beta, pk, rk);
		rkrk = rkp1rkp1;

		//----------Check convergence----------
		T rnorm = sqrt(rkrk);
		//std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}		
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************BiCGSTAB method********************
template<class T>
std::vector<T> BiCGSTAB(CSR<T>& _A, std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Initialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> rdash = rk;
	std::vector<T> pk = rk;
	T rdashrk = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));

	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A*pk;
		T alpha = rdashrk/std::inner_product(rdash.begin(), rdash.end(), Apk.begin(), T());
		std::vector<T> sk = zeaxpby(1.0, rk, -alpha, Apk);
		std::vector<T> Ask = _A*sk;
		T omega = std::inner_product(Ask.begin(), Ask.end(), sk.begin(), T())/std::inner_product(Ask.begin(), Ask.end(), Ask.begin(), T());
		xeaxpbypcz(1.0, xk, alpha, pk, omega, sk);
		rk = zeaxpby(1.0, sk, -omega, Ask);
		T rdashrkp1 = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
		T beta = alpha/omega*rdashrkp1/rdashrk;
		xeaxpbypcz(beta, pk, 1.0, rk, -beta*omega, Apk);
		rdashrk = rdashrkp1;

		//----------Check convergence----------
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		//std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			//std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************BiCGSTAB2 method********************
template<class T>
std::vector<T> BiCGSTAB2(CSR<T>& _A, std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Initialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> rdash = rk;
	std::vector<T> pk(_b.size(), T());
	std::vector<T> uk(_b.size(), T());
	std::vector<T> tkm1(_b.size(), T()); 
	std::vector<T> wk(_b.size(), T());
	std::vector<T> zk(_b.size(), T());
	T beta = T();
	T rdashrk = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));

	//----------Iteration----------
	for(int k = 0; k < _itrmax; k++) {
		xeaxpbypcz(beta, pk, 1.0, rk, -beta, uk);
		std::vector<T> Apk = _A*pk;
		T alpha = rdashrk/std::inner_product(rdash.begin(), rdash.end(), Apk.begin(), T());
		std::vector<T> yk = zeavpbwpcxpdy(1.0, tkm1, -1.0, rk, -alpha, wk, alpha, Apk);
		std::vector<T> tk = zeaxpby(1.0, rk, -alpha, Apk);
		T zeta, ita;
		std::vector<T> Atk = _A*tk;
		T Att = std::inner_product(Atk.begin(), Atk.end(), tk.begin(), T());
		T AtAt = std::inner_product(Atk.begin(), Atk.end(), Atk.begin(), T());
		if(k%2 == 0) {
			zeta = Att/AtAt;
			ita = T();
		} else {
			T yy = std::inner_product(yk.begin(), yk.end(), yk.begin(), T());
			T yt = std::inner_product(yk.begin(), yk.end(), tk.begin(), T());
			T Aty = std::inner_product(Atk.begin(), Atk.end(), yk.begin(), T());
			zeta = (yy*Att - yt*Aty)/(AtAt*yy - Aty*Aty);
			ita = (AtAt*yt - Aty*Att)/(AtAt*yy - Aty*Aty);
		}
		zeawpbxmypcz(zeta, Apk, ita, tkm1, rk, beta, uk);
		xeaxpbypcz(ita, zk, zeta, rk, -alpha, uk);
		xeaxpbypcz(1.0, xk, alpha, pk, 1.0, zk);
		rk = zeawpbxpcy(1.0, tk, -ita, yk, -zeta, Atk);
		T rdashrkp1 = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
		beta = alpha*rdashrkp1/(zeta*rdashrk);
		wk = zeaxpby(1.0, Atk, beta, Apk);
		rdashrk = rdashrkp1;

		tkm1 = tk;
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		//std::cout << "k = " << k << "\teps = " << rnorm/bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			//std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************Incomplete LU(0) decomposition********************
template<class T>
CSR<T> ILU0(CSR<T>& _A) {
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
						} else {
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


//********************Solve with ILU(0)*******************
template<class T>
std::vector<T> PreILU0(CSR<T>& _A, std::vector<T>& _b) {
	//----------Solve Ly=b with ILU(0)----------
	std::vector<T> v = std::vector<T>(_b);
	for (int i = 0; i < _b.size(); i++) {
		for (int k = _A.indptr[i]; k < _A.indptr[i + 1]; k++) {
			if (_A.indices[k] < i) {
				v[i] -= _A.data[k] * v[_A.indices[k]];
			} else {
				break;
			}
		}
	}

	//----------Solve Ux=y with ILU(0)----------
	for (int i = _b.size() - 1; i >= 0; i--) {
		for (int k = _A.indptr[i + 1] - 1; k >= _A.indptr[i]; k--) {
			if (_A.indices[k] > i) {
				v[i] -= _A.data[k] * v[_A.indices[k]];
			} else {
				break;
			}
		}
		v[i] /= _A.get(i, i);
	}

	return v;
}


//*******************ILU(0) preconditioning CG method********************
template<class T>
std::vector<T> ILU0CG(CSR<T>& _A, CSR<T>& _M, const std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Initialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> pk = PreILU0(_M, rk);				//Preconditioning
	std::vector<T> Mrk = pk;							
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));
	T Mrkrk = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());
	
	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A*pk;
		T alpha = Mrkrk/std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
		xexpay(xk, alpha, pk);
		xexpay(rk, -alpha, Apk);
		Mrk = PreILU0(_M, rk);							//Preconditioning
		T Mrkp1rkp1 = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());
		T beta = Mrkp1rkp1/Mrkrk;
		xeaxpy(beta, pk, Mrk);
		Mrkrk = Mrkp1rkp1;

		//----------Check convergence----------
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		//std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//*******************ILU(0) preconditioning BiCGSTAB method*******************
template<class T>
std::vector<T> ILU0BiCGSTAB(CSR<T>& _A, CSR<T>& _M, std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Iniialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> rdash = rk;
	std::vector<T> pk = rk;
	T rdashrk = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));

	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Mpk = PreILU0(_M, pk);		//Preconditioning
		std::vector<T> AMpk = _A*Mpk;
		T alpha = rdashrk/std::inner_product(rdash.begin(), rdash.end(), AMpk.begin(), T());
		std::vector<T> sk = zeaxpby(1.0, rk, -alpha, AMpk);
		std::vector<T> Msk = PreILU0(_M, sk);		//Preconditioning
		std::vector<T> AMsk = _A*Msk;
		T omega = std::inner_product(AMsk.begin(), AMsk.end(), sk.begin(), T())/std::inner_product(AMsk.begin(), AMsk.end(), AMsk.begin(), T());
		xeaxpbypcz(1.0, xk, alpha, Mpk, omega, Msk);
		rk = zeaxpby(1.0, sk, -omega, AMsk);
		T rdashrkp1 = std::inner_product(rdash.begin(), rdash.end(), rk.begin(), T());
		T beta = alpha/omega*rdashrkp1/rdashrk;
		xeaxpbypcz(beta, pk, 1.0, rk, -beta*omega, AMpk);
		rdashrk = rdashrkp1;

		//----------Check convergence----------
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************Get diagonal vector of matrix _A********************
template<class T>
std::vector<T> GetDiagonal(CSR<T>& _A) {
	std::vector<T> v(_A.ROWS);
	for (int i = 0; i < _A.ROWS; i++) {
		v[i] = _A.get(i, i);
	}
	return v;
}


//********************Scaling matrix********************
template<class T>
std::vector<T> Scaling(std::vector<T>& _D, std::vector<T>& _b) {
	std::vector<T> v(_D.size());
	for (int i = 0; i < _D.size(); i++) {
		v[i] = _b[i] / _D[i];
	}
	return v;
}


//********************Scaling preconditioning CG method********************
template<class T>
std::vector<T> ScalingCG(CSR<T>& _A, const std::vector<T>& _b, int _itrmax, T _eps) {
	//----------Initialize----------
	std::vector<T> D = GetDiagonal(_A);				//Scaling A matrix
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> pk = Scaling(D, rk);				//Scaling rk
	std::vector<T> Mrk = pk;
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));
	T Mrkrk = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());

	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A*pk;
		T alpha = Mrkrk/std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
		xexpay(xk, alpha, pk);
		xexpay(rk, -alpha, Apk);
		Mrk = Scaling(D, rk);						//Scaling rkp1
		T Mrkp1rkp1 = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());
		T beta = Mrkp1rkp1/Mrkrk;
		xeaxpy(beta, pk, Mrk);
		Mrkrk = Mrkp1rkp1;

		//----------Check convergence----------
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		//std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}


//********************Solve with SOR********************
template<class T>
std::vector<T> SOR(CSR<T>& _A, std::vector<T>& _b, T _w, int _itrmax, T _eps) {
	std::vector<T> x = std::vector<T>(_A.ROWS, T());
	T error = T();
	for (int itr = 0; itr < _itrmax; itr++) {
		error = T();
		for (int i = 0; i < _A.ROWS; i++) {
			T Aii = T();
			T tmp = x[i];
			x[i] = _b[i];
			for (int k = _A.indptr[i]; k < _A.indptr[i + 1]; k++) {
				int j = _A.indices[k];
				if (i != j) {
					x[i] -= _A.data[k] * x[j];
				} else {
					Aii = _A.data[k];
				}
			}
			x[i] /= Aii;
			x[i] = tmp + _w * (x[i] - tmp);

			error += fabs((tmp - x[i]) / tmp);
		}
		//std::cout << error << std::endl;
		if (error < _eps) {
			//std::cout << "\tConvergence:" << itr << std::endl;
			return x;
		}
	}

	//std::cout << "\t" << error << std::endl;
	return x;
}


//********************SOR preconditioning CG method********************
template<class T>
std::vector<T> SORCG(CSR<T>& _A, const std::vector<T>& _b, int _itrmax, T _eps, T _soromega, int _soritermax, T _soreps) {
	//----------Initialize----------
	std::vector<T> xk(_b.size(), T());
	std::vector<T> rk = subtract(_b, _A*xk);
	std::vector<T> pk = SOR(_A, rk, _soromega, _soritermax, _soreps);		//Preconditioning SOR
	std::vector<T> Mrk = pk;
	T bnorm = sqrt(std::inner_product(_b.begin(), _b.end(), _b.begin(), T()));
	T Mrkrk = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());

	//----------Iteration----------
	for (int k = 0; k < _itrmax; ++k) {
		std::vector<T> Apk = _A*pk;
		T alpha = Mrkrk/std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
		xexpay(xk, alpha, pk);
		xexpay(rk, -alpha, Apk);
		Mrk = SOR(_A, rk, _soromega, _soritermax, _soreps);					//Preconditioning SOR
		T Mrkp1rkp1 = std::inner_product(Mrk.begin(), Mrk.end(), rk.begin(), T());
		T beta = Mrkp1rkp1/Mrkrk;
		xeaxpy(beta, pk, Mrk);
		Mrkrk = Mrkp1rkp1;

		//----------Check convergence----------
		T rnorm = sqrt(std::inner_product(rk.begin(), rk.end(), rk.begin(), T()));
		//std::cout << "k = " << k << "\teps = " << rnorm / bnorm << std::endl;
		if (rnorm < _eps*bnorm) {
			std::cout << "\tConvergence:" << k << std::endl;
			return xk;
		}
	}

	std::cout << "\nConvergence:faild" << std::endl;
	return xk;
}