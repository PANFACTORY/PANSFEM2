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
    std::vector<T> qk = std::vector<T>(n, 1.0/sqrt((T)n));

    for(int k = 0; k < n; k++){
        std::vector<T> v = _A*qk;
        alpha[k] = std::inner_product(qk.begin(), qk.end(), v.begin(), T());
        Subtruct(v, beta[k - 1], qkm1, alpha[k], qk);
        beta[k] = std::inner_product(v.begin(), v.end(), v.begin(), T());
        qkm1 = qk;
        std::transform(v.begin(), v.end(), qk.begin(), [=](T _vi) { return _vi / beta[k]; });
    }
}


