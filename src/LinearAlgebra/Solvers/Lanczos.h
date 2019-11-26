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


//********************Vector subtruct********************
template<class T>
void Subtruct(std::vector<T>& _v, T _beta, const std::vector<T>& _qkm1, T _alpha, const std::vector<T>& _qk){
    auto qkm1i = _qkm1.begin(), qki = _qk.begin();
	for (auto &vi : _v) {
		vi -= _beta*(*qkm1i) + _alpha*(*qki);
		++qkm1i;
		++qki;
	}
} 


//********************Lanczos process********************
template<class T>
void LanczosProcess(CSR<T>& _A, std::vector<T>& _alpha, std::vector<T>& _beta){
    int n = _A.ROWS;                    //Degree of matrix
    
    _alpha = std::vector<T>(n);         //Values of diagonal
    _beta = std::vector<T>(n);          //Values of side of diagonal

    std::vector<T> qkm1 = std::vector<T>(n, T());
    std::vector<T> qk = std::vector<T>(n, T());
    qk[0] = 1.0;

    for(int k = 0; k < n; k++){
        std::vector<T> v = _A*qk;
        _alpha[k] = std::inner_product(qk.begin(), qk.end(), v.begin(), T());
        if(k == 0) {
            Subtruct(v, T(), qkm1, _alpha[k], qk);
        } else {
            Subtruct(v, _beta[k - 1], qkm1, _alpha[k], qk);
        }
        _beta[k] = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), T()));
        qkm1 = qk;
        std::transform(v.begin(), v.end(), qk.begin(), [=](T _vi) { return _vi/_beta[k]; });
    }
}


//********************Bisection method********************
template<class T>
T BisectionMethod(const std::vector<T>& _alpha, const std::vector<T>& _beta, int _mode){
    int n = _alpha.size();
    
    T b = fabs(_alpha[0]) + fabs(_beta[0]);
    for(int i = 1; i < n; i++){
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
        for(int k = 0; k < n - 1; k++){
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
    int n = _alpha.size();
    std::vector<T> yk = std::vector<T>(n, T());
    yk[0] = 1.0;

    for(int k = 0; k < 10; k++){
        //.....Solve [T - lambda*E]{ykp1}={yk} with Thomas method.....
        std::vector<T> p = std::vector<T>(n);
        std::vector<T> q = std::vector<T>(n);
        
        p[0] = -_beta[0]/(_alpha[0] - _lambda);
        q[0] = yk[0]/(_alpha[0] - _lambda);
        for(int i = 1; i < n; i++){
            p[i] = -_beta[i]/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
            q[i] = (yk[i] - _beta[i - 1]*q[i - 1])/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
        }

        yk[n - 1] = q[n - 1];
        for(int i = n - 2; i >= 0; i--){
            yk[i] = p[i]*yk[i + 1] + q[i];
        }

        //.....Normalize yk.....
        T ykNorm = sqrt(std::inner_product(yk.begin(), yk.end(), yk.begin(), T()));
        for(int i = 0; i < n; i++){
            yk[i] /= ykNorm;
        }
    }
    
    return yk;
}