//*****************************************************************************
//  Title       :   src/LinearAlgebra/Solvers/CG.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/11/09
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <iostream>
#include <cassert>

namespace PANSFEM2 {
    //  {x} = a{y}
    template<class T>
    inline void xeay(T *_x, T _a, T *_y, int _n) {
        for (int i = 0; i < _n; i++) {
            _x[i] = _a*_y[i];
        }
    }

    //  {x} = a{y} + b{z}
    template<class T>
    inline void xeaypbz(T *_x, T _a, T *_y, T _b, T *_z, int _n) {
        for (int i = 0; i < _n; ++i) {
            _x[i] = _a*_y[i] + _b*_z[i];
        }
    }

    //  {w} = a{x} + b{y} + c{z}
    template<class T>
    inline void weaxpbypcz(T *_w, T _a, T *_x, T _b, T *_y, T _c, T *_z, int _n) {
        for (int i = 0; i < _n; ++i) {
            _w[i] = _a*_x[i] + _b*_y[i] + _c*_z[i];
        }
    }

    //  {v} = a{w} + b{x} + c{y} + d{z}
    template<class T>
    inline void veawpbxpcypdz(T *_v, T _a, T *_w, T _b, T *_x, T _c, T *_y, T _d, T *_z, int _n) {
        for (int i = 0; i < _n; ++i) {
            _v[i] = _a*_w[i] + _b*_x[i] + _c*_y[i] + _d*_z[i];
        }
    }

    //  Inner product of {x} and {y}
    template<class T>
    inline T innerproduct(T *_x, T *_y, int _n) {
        T innerproduct = T();
        for (int i = 0; i < _n; ++i) {
            innerproduct += _x[i]*_y[i];
        }
        return innerproduct;
    }

    //********************CG solver for dense matrix********************
    template<class T>
    void CG(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *p = new T[_n], *Ap = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            p[i] = r[i];
        }
        T rkrk = innerproduct(r, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);

        //----------Iteration----------
        int k = 0;
        for (; k < _itrmax; ++k) {
            for (int i = 0; i < _n; ++i) {
                Ap[i] = T();
                for (int j = 0; j < _n; ++j) {
                    Ap[i] += _A[_n*i + j]*p[j];
                }
            }   //  [A]{p}
            T alpha = rkrk/innerproduct(p, Ap, _n);
            xeaypbz(_x, 1.0, _x, alpha, p, _n);
            xeaypbz(r, 1.0, r, -alpha, Ap, _n);
            T rkp1rkp1 = innerproduct(r, r, _n);
            T beta = rkp1rkp1/rkrk;
            xeaypbz(p, beta, p, 1.0, r, _n);
            rkrk = rkp1rkp1;

            //----------Check convergence----------
            if (_isdebug) {
                std::cout << "k = " << k << "\teps = " << sqrt(rkrk/bb) << std::endl;    
            }
            if (rkrk < epsepsbb) {
                break;
            }		
        }

        //----------Show result----------
        if (_isdebug) {
            if (k < _itrmax) {
                std::cout << "\tConvergence:" << k << std::endl;
            } else {
                std::cout << "\nConvergence:faild" << std::endl;
            }
        }
        
        delete[] r;
        delete[] p;
        delete[] Ap;
    }

    //********************BiCGSTAB solver for dense matrix********************
    template<class T>
    void BiCGSTAB(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *rdash = new T[_n], *p = new T[_n], *Ap = new T[_n], *s = new T[_n], *As = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            rdash[i] = r[i];
            p[i] = r[i];
        }
        T rdashrk = innerproduct(rdash, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);

        //----------Iteration----------
        int k = 0;
        for (; k < _itrmax; ++k) {
            for (int i = 0; i < _n; ++i) {
                Ap[i] = T();
                for (int j = 0; j < _n; ++j) {
                    Ap[i] += _A[_n*i + j]*p[j];
                }
            }   //  [A]{p}
            T alpha = rdashrk/innerproduct(rdash, Ap, _n);
            xeaypbz(s, 1.0, r, -alpha, Ap, _n);
            for (int i = 0; i < _n; ++i) {
                As[i] = T();
                for (int j = 0; j < _n; ++j) {
                    As[i] += _A[_n*i + j]*s[j];
                }
            }   //  [A]{s}
            T omega = innerproduct(As, s, _n)/innerproduct(As, As, _n);
            weaxpbypcz(_x, 1.0, _x, alpha, p, omega, s, _n);
            xeaypbz(r, 1.0, s, -omega, As, _n);
            T rdashrkp1 = innerproduct(rdash, r, _n);
            T beta = alpha/omega*rdashrkp1/rdashrk;
            weaxpbypcz(p, beta, p, 1.0, r, -beta*omega, Ap, _n);
            rdashrk = rdashrkp1;

            //----------Check convergence----------
            T rkrk = innerproduct(r, r, _n);
            if (_isdebug) {
                std::cout << "k = " << k << "\teps = " << sqrt(rkrk/bb) << std::endl;    
            }
            if (rkrk < epsepsbb) {
                break;
            }	
        }

        //----------Show result----------
        if (_isdebug) {
            if (k < _itrmax) {
                std::cout << "\tConvergence:" << k << std::endl;
            } else {
                std::cout << "\nConvergence:faild" << std::endl;
            }
        }

        delete[] r;
        delete[] rdash;
        delete[] p;
        delete[] Ap;
        delete[] s;
        delete[] As;
    }

    //********************BiCGSTAB2 solver for dense matrix********************
    template<class T>
    void BiCGSTAB2(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *rdash = new T[_n], *p = new T[_n], *Ap = new T[_n], *t = new T[_n], *At = new T[_n], *u = new T[_n], *w = new T[_n], *y = new T[_n], *z = new T[_n], *tkm1 = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            rdash[i] = r[i];
            p[i] = T();
            u[i] = T();
            w[i] = T();
            z[i] = T();
            tkm1[i] = T();
        }
        T beta = T();
        T rdashrk = innerproduct(rdash, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);
        
        //----------Iteration----------
	    int k = 0;
        for(; k < _itrmax; k++) {
            weaxpbypcz(p, beta, p, 1.0, r, -beta, u, _n);
            for (int i = 0; i < _n; ++i) {
                Ap[i] = T();
                for (int j = 0; j < _n; ++j) {
                    Ap[i] += _A[_n*i + j]*p[j];
                }
            }   //  [A]{p}
            T alpha = rdashrk/innerproduct(rdash, Ap, _n);
            veawpbxpcypdz(y, 1.0, tkm1, -1.0, r, -alpha, w, alpha, Ap, _n);
            xeaypbz(t, 1.0, r, -alpha, Ap, _n);
            T zeta, ita;
            for (int i = 0; i < _n; ++i) {
                At[i] = T();
                for (int j = 0; j < _n; ++j) {
                    At[i] += _A[_n*i + j]*t[j];
                }
            }   //  [A]{t}
            T Att = innerproduct(At, t, _n);
            T AtAt = innerproduct(At, At, _n);
            if (k%2 == 0) {
                zeta = Att/AtAt;
                ita = T();
            } else {
                T yy = innerproduct(y, y, _n);
                T yt = innerproduct(y, t, _n);
                T Aty = innerproduct(At, y, _n);
                zeta = (yy*Att - yt*Aty)/(AtAt*yy - Aty*Aty);
			    ita = (AtAt*yt - Aty*Att)/(AtAt*yy - Aty*Aty);
            }
            veawpbxpcypdz(u, zeta, Ap, ita, tkm1, -ita, r, beta*ita, u, _n);
            weaxpbypcz(z, ita, z, zeta, r, -alpha, u, _n);
            weaxpbypcz(_x, 1.0, _x, alpha, p, 1.0, z, _n);
            weaxpbypcz(r, 1.0, t, -ita, y, -zeta, At, _n);
            T rdashrkp1 = innerproduct(rdash, r, _n);
            beta = alpha*rdashrkp1/(zeta*rdashrk);
            xeaypbz(w, 1.0, At, beta, Ap, _n);
            rdashrk = rdashrkp1;
            xeay(tkm1, 1.0, t, _n);

            //----------Check convergence----------
            T rkrk = innerproduct(r, r, _n);
            if (_isdebug) {
                std::cout << "k = " << k << "\teps = " << sqrt(rkrk/bb) << std::endl;    
            }
            if (rkrk < epsepsbb) {
                break;
            }
        }    

        //----------Show result----------
        if (_isdebug) {
            if (k < _itrmax) {
                std::cout << "\tConvergence:" << k << std::endl;
            } else {
                std::cout << "\nConvergence:faild" << std::endl;
            }
        }

        delete[] r;
        delete[] rdash;
        delete[] p;
        delete[] Ap;
        delete[] t;
        delete[] At;
        delete[] u;
        delete[] w;
        delete[] y;
        delete[] z;
        delete[] tkm1;
    }
}