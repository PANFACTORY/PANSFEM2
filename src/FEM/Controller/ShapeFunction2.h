//*****************************************************************************
//Title		:src/FEM/Controller/ShapeFunction2.h
//Author	:Tanabe Yuta
//Date		:2020/05/18
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace PANSFEM2 {
    //********************Lagrange interpolation********************
    template<class T>
    std::vector<T> LagrangeInterpolation(std::vector<T> _xs, T _x) {
        std::vector<T> N = std::vector<T>(_xs.size(), 1.0);
        for(int i = 0; i < _xs.size(); i++) {
            for(int j = 0; j < _xs.size(); j++) {
                if(i != j) {
                    N[i] *= (_x - _xs[j])/(T)(_xs[i] - _xs[j]);
                }
            }
        }
        return N;
    }


    template<class T>
    std::vector<T> LagrangeInterpolationDerivative(std::vector<T> _xs, T _x) {
        std::vector<T> dNdr = std::vector<T>(_xs.size(), T());
        for(int i = 0; i < _xs.size(); i++) {
            for(int j = 0; j < _xs.size(); j++) {
                if(j != i) {
                    T p = 1.0;
                    for(int k = 0; k < _xs.size(); k++) {
                        if(k != j && k != i) {
                            p *= (_x - _xs[k])/(T)(_xs[i] - _xs[k]);
                        } else if(k != i) {
                            p *= 1.0/(T)(_xs[i] - _xs[k]); 
                        }
                    }
                    dNdr[i] += p;
                }
            }
        }
        return dNdr;
    }
}