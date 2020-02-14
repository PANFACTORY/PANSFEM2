//*****************************************************************************
//  Title		:   LinearAlgebra/Solvers/LU.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/10/02
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <numeric>
#include <cassert>


#include "../Models/Matrix.h"
#include "../Models/Vector.h"


namespace PANSFEM2{
    //**********LU decomposition with Pivotting**********
    template<class T>
    void LU(Matrix<T>& _A, std::vector<int>& _pivot){
        int n = _A.ROW();
        assert(n == _A.COL());
        assert(n == _pivot.size());

        std::iota(_pivot.begin(), _pivot.end(), 0);

        for(int k = 0; k < n - 1; k++){
            //----------Get pivot----------
            T pivv = fabs(_A(k, k));
            int pivi = k;
            for(int i = k + 1; i < n; i++){
                if(pivv < fabs(_A(i, k))){
                    pivv = fabs(_A(i, k));
                    pivi = i;
                }
            }

            //----------Exchange row----------
            if(k != pivi){
                std::swap(_pivot[k], _pivot[pivi]);
                for(int j = 0; j < n; j++){
                    std::swap(_A(k, j), _A(pivi, j));
                }
            }
            
            //----------LU decomposition----------
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


    //**********Solve LU**********
    template<class T>
    void SolveLU(Matrix<T>& _LU, Vector<T>& _b, std::vector<int>& _pivot){
        int n = _LU.ROW();
        assert(n == _LU.COL());
        assert(n == _b.SIZE());
        assert(n == _pivot.size());

        //----------Exchange row of b----------
        Vector<T> tmp = _b;
        for(int i = 0; i < n; i++){
            _b(i) = tmp(_pivot[i]);
        }
        
        //----------Solve Ly=b----------
        for(int i = 0; i < n; i++){
            for(int j = 0; j < i; j++){
                _b(i) -= _LU(i, j)*_b(j);
            }
        }

        //----------Solve Ux=y----------
        for(int i = n - 1; i >= 0; i--){
            for(int j = i + 1; j < n; j++){
                _b(i) -= _LU(i, j)*_b(j);
            }
            _b(i) /= _LU(i, i);
        }
    }
}