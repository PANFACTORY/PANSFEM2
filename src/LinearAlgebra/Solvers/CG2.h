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
        T *r = new T[_n], *rdash = new T[_n], *p = new T[_n], *Ap = new T[_n], *s = new T[_n], *As = new T[_n], *u = new T[_n], *w = new T[_n], *y = new T[_n], *z = new T[_n], *skm1 = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            rdash[i] = r[i];
            p[i] = r[i];
            u[i] = T();
            w[i] = T();
            z[i] = T();
            skm1[i] = T();
        }
        T beta = T();
        T rdashrk = innerproduct(rdash, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);
        
        //----------Iteration----------
	    int k = 0;
        for(; k < _itrmax; k++) {
            for (int i = 0; i < _n; ++i) {
                Ap[i] = T();
                for (int j = 0; j < _n; ++j) {
                    Ap[i] += _A[_n*i + j]*p[j];
                }
            }   //  [A]{p}
            T alpha = rdashrk/innerproduct(rdash, Ap, _n);
            veawpbxpcypdz(y, 1.0, skm1, -1.0, r, -alpha, w, alpha, Ap, _n);
            xeaypbz(s, 1.0, r, -alpha, Ap, _n);
            for (int i = 0; i < _n; ++i) {
                As[i] = T();
                for (int j = 0; j < _n; ++j) {
                    As[i] += _A[_n*i + j]*s[j];
                }
            }   //  [A]{s}
            T Ass = innerproduct(As, s, _n);
            T AsAs = innerproduct(As, As, _n);
            T omega = Ass/AsAs, ita = T();
            if (k%2 != 0) {
                T yy = innerproduct(y, y, _n);
                T ys = innerproduct(y, s, _n);
                T Asy = innerproduct(As, y, _n);
                omega = (yy*Ass - ys*Asy)/(AsAs*yy - Asy*Asy);
			    ita = (AsAs*ys - Asy*Ass)/(AsAs*yy - Asy*Asy);
            }
            veawpbxpcypdz(u, omega, Ap, ita, skm1, -ita, r, beta*ita, u, _n);
            weaxpbypcz(z, ita, z, omega, r, -alpha, u, _n);
            weaxpbypcz(_x, 1.0, _x, alpha, p, 1.0, z, _n);
            weaxpbypcz(r, 1.0, s, -omega, As, -ita, y, _n);
            T rdashrkp1 = innerproduct(rdash, r, _n);
            beta = alpha/omega*rdashrkp1/rdashrk;
            xeaypbz(w, 1.0, As, beta, Ap, _n);
            weaxpbypcz(p, beta, p, 1.0, r, -beta, u, _n);
            rdashrk = rdashrkp1;
            xeay(skm1, 1.0, s, _n);

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
        delete[] u;
        delete[] w;
        delete[] y;
        delete[] z;
        delete[] skm1;
    }

    //********************GPBi-CG solver for dense matrix********************
    template<class T>
    void GPBiCG(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20) {
        
    }

    //  Solve [D]{x} = {b}
    template<class T>
    inline void SolveDxey(T *_D, T *_x, T *_y, int _n, T _eps = 1.0e-20) {
        assert(0 < _n && T() < _eps);
        for (int i = 0; i < _n; i++) {
            _x[i] = _eps < fabs(_D[i]) ? _y[i]/_D[i] : _y[i];
        }
    }

    //********************Scaling preconditioned CG solver for dense matrix********************
    template<class T>
    void ScalingCG(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20, T _scalingmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *p = new T[_n], *Ap = new T[_n], *D = new T[_n], *Mr = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            D[i] = _A[_n*i + i];
            Mr[i] = _scalingmin < fabs(D[i]) ? r[i]/D[i] : r[i];
            p[i] = Mr[i];
        }
        T Mrkrk = innerproduct(Mr, r, _n);
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
            T alpha = Mrkrk/innerproduct(p, Ap, _n);
            xeaypbz(_x, 1.0, _x, alpha, p, _n);
            xeaypbz(r, 1.0, r, -alpha, Ap, _n);
            SolveDxey(D, Mr, r, _n, _scalingmin);
            T Mrkp1rkp1 = innerproduct(Mr, r, _n);
            T beta = Mrkp1rkp1/Mrkrk;
            xeaypbz(p, beta, p, 1.0, Mr, _n);
            Mrkrk = Mrkp1rkp1;

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
        delete[] p;
        delete[] Ap;
        delete[] D;
        delete[] Mr;
    }

    //********************Scaling preconditioned BiCGSTAB solver for dense matrix********************
    template<class T>
    void ScalingBiCGSTAB(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20, T _scalingmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *rdash = new T[_n], *p = new T[_n], *Mp = new T[_n], *AMp = new T[_n], *s = new T[_n], *Ms = new T[_n], *AMs = new T[_n], *D = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            D[i] = _A[_n*i + i];
            rdash[i] = r[i];
            p[i] = r[i];
        }
        T rdashrk = innerproduct(rdash, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);

        //----------Iteration----------
        int k = 0;
        for (; k < _itrmax; ++k) {
            SolveDxey(D, Mp, p, _n, _scalingmin);
            for (int i = 0; i < _n; ++i) {
                AMp[i] = T();
                for (int j = 0; j < _n; ++j) {
                    AMp[i] += _A[_n*i + j]*Mp[j];
                }
            }   //  [A][M]{p}
            T alpha = rdashrk/innerproduct(rdash, AMp, _n);
            xeaypbz(s, 1.0, r, -alpha, AMp, _n);
            SolveDxey(D, Ms, s, _n, _scalingmin);
            for (int i = 0; i < _n; ++i) {
                AMs[i] = T();
                for (int j = 0; j < _n; ++j) {
                    AMs[i] += _A[_n*i + j]*Ms[j];
                }
            }   //  [A][M]{s}
            T omega = innerproduct(AMs, s, _n)/innerproduct(AMs, AMs, _n);
            weaxpbypcz(_x, 1.0, _x, alpha, Mp, omega, Ms, _n);
            xeaypbz(r, 1.0, s, -omega, AMs, _n);
            T rdashrkp1 = innerproduct(rdash, r, _n);
            T beta = alpha/omega*rdashrkp1/rdashrk;
            weaxpbypcz(p, beta, p, 1.0, r, -beta*omega, AMp, _n);
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
        delete[] Mp;
        delete[] AMp;
        delete[] s;
        delete[] Ms;
        delete[] AMs;
        delete[] D;
    }

    //********************Scaling preconditioned BiCGSTAB2 solver for dense matrix********************
    template<class T>
    void ScalingBiCGSTAB2(T *_A, T *_b, T *_x, int _n, int _itrmax, T _eps, bool _isdebug = false, T _bbmin = 1.0e-20, T _scalingmin = 1.0e-20) {
        assert(0 < _n && 0 < _itrmax && T() < _eps);

        //----------Initialize----------
        T *r = new T[_n], *rdash = new T[_n], *p = new T[_n], *Mp = new T[_n], *AMp = new T[_n], *s = new T[_n], *Ms = new T[_n], *AMs = new T[_n], *u = new T[_n], *w = new T[_n], *y = new T[_n], *z = new T[_n], *Mz = new T[_n], *skm1 = new T[_n], *D = new T[_n];
        for (int i = 0; i < _n; ++i) {
            r[i] = _b[i];
            for (int j = 0; j < _n; ++j) {
                r[i] -= _A[_n*i + j]*_x[j];
            }
            D[i] = _A[_n*i + i];
            rdash[i] = r[i];
            p[i] = r[i];
            u[i] = T();
            w[i] = T();
            z[i] = T();
            skm1[i] = T();
        }
        T beta = T();
        T rdashrk = innerproduct(rdash, r, _n);
        T bb = innerproduct(_b, _b, _n);
        T epsepsbb = _eps*_eps*std::max(bb, _bbmin);
        
        //----------Iteration----------
	    int k = 0;
        for(; k < _itrmax; k++) {
            SolveDxey(D, Mp, p, _n, _scalingmin);
            for (int i = 0; i < _n; ++i) {
                AMp[i] = T();
                for (int j = 0; j < _n; ++j) {
                    AMp[i] += _A[_n*i + j]*Mp[j];
                }
            }   //  [A][M]{p}
            T alpha = rdashrk/innerproduct(rdash, AMp, _n);
            veawpbxpcypdz(y, 1.0, skm1, -1.0, r, -alpha, w, alpha, AMp, _n);
            xeaypbz(s, 1.0, r, -alpha, AMp, _n);
            SolveDxey(D, Ms, s, _n, _scalingmin);
            for (int i = 0; i < _n; ++i) {
                AMs[i] = T();
                for (int j = 0; j < _n; ++j) {
                    AMs[i] += _A[_n*i + j]*Ms[j];
                }
            }   //  [A][M]{s}
            T AMss = innerproduct(AMs, s, _n);
            T AMsAMs = innerproduct(AMs, AMs, _n);
            T omega = AMss/AMsAMs, ita = T();
            if (k%2 != 0) {
                T yy = innerproduct(y, y, _n);
                T ys = innerproduct(y, s, _n);
                T AMsy = innerproduct(AMs, y, _n);
                omega = (yy*AMss - ys*AMsy)/(AMsAMs*yy - AMsy*AMsy);
			    ita = (AMsAMs*ys - AMsy*AMss)/(AMsAMs*yy - AMsy*AMsy);
            }
            veawpbxpcypdz(u, omega, AMp, ita, skm1, -ita, r, beta*ita, u, _n);
            weaxpbypcz(z, ita, z, omega, r, -alpha, u, _n);
            SolveDxey(D, Mz, z, _n, _scalingmin);
            weaxpbypcz(_x, 1.0, _x, alpha, Mp, 1.0, Mz, _n);
            weaxpbypcz(r, 1.0, s, -omega, AMs, -ita, y, _n);
            T rdashrkp1 = innerproduct(rdash, r, _n);
            beta = alpha/omega*rdashrkp1/rdashrk;
            xeaypbz(w, 1.0, AMs, beta, AMp, _n);
            weaxpbypcz(p, beta, p, 1.0, r, -beta, u, _n);
            rdashrk = rdashrkp1;
            xeay(skm1, 1.0, s, _n);

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
        delete[] Mp;
        delete[] AMp;
        delete[] s;
        delete[] Ms;
        delete[] AMs;
        delete[] u;
        delete[] w;
        delete[] y;
        delete[] z;
        delete[] Mz;
        delete[] skm1;
        delete[] D;
    }
}