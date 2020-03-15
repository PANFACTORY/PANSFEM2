//*****************************************************************************
//Title		:src/FEM/Equation/Numeric.h
//Author	:Tanabe Yuta
//Date		:2020/03/15
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Numerical integration on line********************
	template<class T, template<class>class SF, template<class>class IC, class F>
	T NumericIntegrationOnLine(std::vector<Vector<T> >& _x, std::vector<int>& _element, F _f) {
        T value = T();

        Matrix<T> X = Matrix<T>(0, 1);
		for(auto i : _element){
			X = X.Vstack(_x[i].Transpose());
		}

        for(int g = 0; g < IC<T>::N; g++){
            //----------Get shapefunction----------
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get x and f----------
			Vector<T> x = X.Transpose()*N;
            Matrix<T> dXdr = dNdr*X;
			T dl = sqrt((dXdr*dXdr.Transpose())(0, 0));

            //----------Update value----------
            value += IC<T>::Weights[g][0]*_f(x)*dl;
        }

        return value;
    }
}