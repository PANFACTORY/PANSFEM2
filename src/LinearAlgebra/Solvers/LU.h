//*****************************************************************************
//  Title		:   LinearAlgebra/Solvers/LU.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/10/02
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


#include "../Models/Matrix.h"


namespace PANSFEM2{
    //**********LU decomposition**********
    template<class T>
    void LU(Matrix<T>& _A){
        int n = _A.ROW();
        assert(n == _A.COL());

        for(int k = 0; k < n - 1; k++){
            for(int i = k + 1; i < n; i++){
                _A(i, k) /= _A(k, k); 
            }
            for(int j = k + 1; j < n; j++){
                for(int i = k + 1; i < n; i++){
                    _A(i, j) -= _A(i, k)*_A(k, j);
                }
            }
        }
    }
}