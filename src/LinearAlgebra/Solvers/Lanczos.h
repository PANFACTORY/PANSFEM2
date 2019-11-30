//*****************************************************************************
//  Title       :src/LinearAlgebra/Solvers/Lanczos.h
//  Author      :Tanabe Yuta
//  Date        :2019/11/22
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <numeric>


#include "../Models/CSR.h"
#include "CG.h"


//********************Lanczos process********************
template<class T>
void LanczosProcess(CSR<T>& _A, std::vector<T>& _alpha, std::vector<T>& _beta, std::vector<std::vector<T> >& _q, int _m){
    _alpha = std::vector<T>(_m);        //Values of diagonal
    _beta = std::vector<T>(_m);         //Values of side of diagonal
    _q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));      //Orthogonal vectors
    for(int i = 0; i < _A.ROWS; i++){
        _q[0][i] = _A.get(i, i);
    }
    T q0Norm = sqrt(std::inner_product(_q[0].begin(), _q[0].end(), _q[0].begin(), T()));
    std::transform(_q[0].begin(), _q[0].end(), _q[0].begin(), [=](T _q0i) { return _q0i/q0Norm; });

    for(int k = 0; k < _m; k++){
        std::vector<T> p = _A*_q[k];
        _alpha[k] = std::inner_product(_q[k].begin(), _q[k].end(), p.begin(), T());
        std::transform(p.begin(), p.end(), _q[k].begin(), p.begin(), [=](T _pi, T _qki) {return _pi - _alpha[k]*_qki; });
        if( k != 0){
            std::transform(p.begin(), p.end(), _q[k - 1].begin(), p.begin(), [=](T _pi, T _qkm1i) {return _pi - _beta[k - 1]*_qkm1i; });
        }
        _beta[k] = sqrt(std::inner_product(p.begin(), p.end(), p.begin(), T()));
        if(k != _m - 1){
            std::transform(p.begin(), p.end(), _q[k + 1].begin(), [=](T _pi) { return _pi/_beta[k]; });
        }
    }
}


//********************Lanczos Inverse Power process********************
template<class T>
void LanczosInversePowerProcess(CSR<T>& _A, std::vector<T>& _alpha, std::vector<T>& _beta, std::vector<std::vector<T> >& _q, int _m){
    int n = _A.ROWS;
    _alpha = std::vector<T>(_m);        //Values of diagonal
    _beta = std::vector<T>(_m);         //Values of side of diagonal
    _q = std::vector<std::vector<T> >(_m, std::vector<T>(n, T()));      //Orthogonal vectors
    for(int i = 0; i < n; i++){
        _q[0][i] = _A.get(i, i);
    }
    T q0Norm = sqrt(std::inner_product(_q[0].begin(), _q[0].end(), _q[0].begin(), T()));
    std::transform(_q[0].begin(), _q[0].end(), _q[0].begin(), [=](T _q0i) { return _q0i/q0Norm; });
    
    int itrmax = std::max(n, 1000);
    for(int k = 0; k < _m; k++){
        std::vector<T> p = ScalingCG(_A, _q[k], itrmax, 1.0e-10);
        if(k != 0){
            std::transform(p.begin(), p.end(), _q[k - 1].begin(), p.begin(), [=](T _pi, T _qkm1i) {return _pi - _beta[k - 1]*_qkm1i; });
        }
        _alpha[k] = std::inner_product(p.begin(), p.end(), _q[k].begin(), T());
        std::transform(p.begin(), p.end(), _q[k].begin(), p.begin(), [=](T _pi, T _qki) {return _pi - _alpha[k]*_qki; });
        _beta[k] = sqrt(std::inner_product(p.begin(), p.end(), p.begin(), T()));
        if(k != _m - 1){
            std::transform(p.begin(), p.end(), _q[k + 1].begin(), [=](T _pi) { return _pi/_beta[k]; });
        }       
    }
}


//********************Lanczos Inverse Power process for General eigenvalue problem********************
template<class T>
void LanczosInversePowerProcessForGeneral(CSR<T>& _A, CSR<T>& _B, std::vector<T>& _alpha, std::vector<T>& _beta, std::vector<std::vector<T> >& _q, int _m){
    int n = _A.ROWS;
    _alpha = std::vector<T>(_m);        //Values of diagonal
    _beta = std::vector<T>(_m);         //Values of side of diagonal
    _q = std::vector<std::vector<T> >(_m, std::vector<T>(n, T()));      //Orthogonal vectors
    for(int i = 0; i < n; i++){
        _q[0][i] = _A.get(i, i);
    }
    T q0Norm = sqrt(std::inner_product(_q[0].begin(), _q[0].end(), _q[0].begin(), T()));
    std::transform(_q[0].begin(), _q[0].end(), _q[0].begin(), [=](T _q0i) { return _q0i/q0Norm; });

    std::vector<T> p = _B*_q[0];
    int itrmax = std::max(n, 1000);
    for(int k = 0; k < _m; k++){
        std::vector<T> s = ScalingCG(_A, p, itrmax, 1.0e-10);
        if(k != 0){
            std::transform(s.begin(), s.end(), _q[k - 1].begin(), s.begin(), [=](T _si, T _qkm1i) {return _si - _beta[k - 1]*_qkm1i; });
        }
        _alpha[k] = std::inner_product(p.begin(), p.end(), s.begin(), T());
        std::transform(s.begin(), s.end(), _q[k].begin(), s.begin(), [=](T _si, T _qki) {return _si - _alpha[k]*_qki; });
        std::vector<T> r = _B*s;
        _beta[k] = sqrt(std::inner_product(r.begin(), r.end(), s.begin(), T()));
        std::transform(r.begin(), r.end(), p.begin(), [=](T _ri) { return _ri/_beta[k]; });
        if(k != _m - 1){
            std::transform(s.begin(), s.end(), _q[k + 1].begin(), [=](T _si) { return _si/_beta[k]; });
        }       
    }
}


//********************Bisection method********************
template<class T>
T BisectionMethod(const std::vector<T>& _alpha, const std::vector<T>& _beta, int _mode){
    int m = _alpha.size();
    
    T b = fabs(_alpha[0]) + fabs(_beta[0]);
    for(int i = 1; i < m; i++){
        if(b < fabs(_beta[i - 1]) + fabs(_alpha[i]) + fabs(_beta[i])){
            b = fabs(_beta[i - 1]) + fabs(_alpha[i]) + fabs(_beta[i]);
        }
    }

    T a1 = -b, a2 = b, lambda;
    while((a2 - a1) > 1.0e-13*b){
        //----------Update lambda----------
        lambda = 0.5*(a1 + a2);

        //----------Count sign change----------
        T pkm1 = 1.0, pk = _alpha[0] - lambda, sign = 1.0;
        int Nlambda = 0;
        if(sign*pk < T()){
            Nlambda++;
            sign *= -1.0;
        }
        for(int k = 0; k < m - 1; k++){
            T pkp1 = (_alpha[k + 1] - lambda)*pk - pow(_beta[k], 2.0)*pkm1;
            if(sign*pkp1 < T()){
                Nlambda++;
                sign *= -1.0;
            }
            pkm1 = pk;
            pk = pkp1;
        }

        //----------Update a1 and a2----------
        if(Nlambda <= _mode){
            a1 = lambda;
        } else {
            a2 = lambda;
        }
    }

    return lambda;
}


//********************Inverse Power method********************
template<class T>
std::vector<T> InversePowerMethod(const std::vector<T>& _alpha, const std::vector<T>& _beta, T _lambda){
    int m = _alpha.size();
    std::vector<T> yk = std::vector<T>(m, T());
    yk[0] = 1.0;

    std::vector<T> p = std::vector<T>(m);
    std::vector<T> q = std::vector<T>(m);
    for(int k = 0; k < 10000; k++){
        //.....Solve [T - lambda*E]{ykp1}={yk} with Thomas method.....
        p[0] = -_beta[0]/(_alpha[0] - _lambda);
        q[0] = yk[0]/(_alpha[0] - _lambda);
        for(int i = 1; i < m; i++){
            p[i] = -_beta[i]/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
            q[i] = (yk[i] - _beta[i - 1]*q[i - 1])/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
        }

        yk[m - 1] = q[m - 1];
        for(int i = m - 2; i >= 0; i--){
            yk[i] = p[i]*yk[i + 1] + q[i];
        }

        //.....Normalize yk.....
        T ykNorm = sqrt(std::inner_product(yk.begin(), yk.end(), yk.begin(), T()));
        std::transform(yk.begin(), yk.end(), yk.begin(), [=](T _yki) { return _yki/ykNorm; });
    }
    
    return yk;
}


//********************Reconvert vector********************
template<class T>
std::vector<T> ReconvertVector(const std::vector<T>& _y, std::vector<std::vector<T> >& _q){
    int n = _q[0].size();
    int m = _y.size();
    std::vector<T> x = std::vector<T>(n, T());

    for(int k = 0; k < m; k++){
        for(int i = 0; i < n; i++){
            x[i] += _q[k][i]*_y[k];
        } 
    }

    return x;
}