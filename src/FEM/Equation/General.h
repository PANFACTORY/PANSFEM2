//*****************************************************************************
//  Title		:   src/FEM/Equation/General.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/04/10
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //**********Get element's area (Plane integration in region)**********
    template<class T, template<class>class SF, template<class>class IC>
    T Area(std::vector<Vector<T> >& _x, std::vector<int>& _element) {
        T area = T();

        //----------Generate cordinate matrix X----------
		Matrix<T> X = Matrix<T>(0, 2);
		for(int i = 0; i < _element.size(); i++){
			X = X.Vstack(_x[_element[i]].Transpose());
		}

        //----------Loop of Gauss Integration----------
        for (int g = 0; g < IC<T>::N; g++) {
			//----------Get shape function and difference of shape function----------
			Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);

			//----------Get difference of shape function----------
			Matrix<T> dXdr = dNdr*X;
			area += dXdr.Determinant()*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }

        return area;
    }


    //**********Get element's center of gravity**********
    template<class T>
    Vector<T> CenterOfGravity(std::vector<Vector<T> >& _x, std::vector<int>& _element) {
        Vector<double> center = Vector<double>(2);
        for(auto i : _element){
            center += _x[i];
        }
        return center/(double)_element.size();
    }


    //**********Get element vector**********
    template<class T>
    Vector<T> ElementVector(std::vector<Vector<T> >& _u, const std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element) {
        int vectorsize = 0;
        for(auto nodetoelementi : _nodetoelement) {
            vectorsize += nodetoelementi.size();
        }

        Vector<T> elementvector = Vector<T>(vectorsize);
        for(int i = 0; i < _nodetoelement.size(); i++) {
            for(auto dou : _nodetoelement[i]) {
                elementvector(dou.second) = _u[_element[i]](dou.first);
            }
        }

        return elementvector;
    }


    //**********Get element vector**********
    template<class T>
    Vector<T> ElementVector(std::vector<Vector<T> >& _u, const std::vector<std::vector<std::vector<std::pair<int, int> > > >& _nodetoelements, const std::vector<std::vector<int> >& _elements) {
        int vectorsize = 0;
        for(auto nodetoelement : _nodetoelements) {
            for(auto nodetoelementi : nodetoelement) {
                vectorsize += nodetoelementi.size();
            }
        }
        
        Vector<T> elementvector = Vector<T>(vectorsize);
        for(int i = 0; i < _nodetoelements.size(); i++) {
            for(int j = 0; j < _nodetoelements[i].size(); j++) {
                for(auto dou : _nodetoelements[i][j]) {
                    elementvector(dou.second) = _u[_elements[i][j]](dou.first);
                }
            }
        }
        
        return elementvector;
    }
}