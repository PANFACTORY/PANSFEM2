//*****************************************************************************
//  Title		:   src/FEM/Equation/Truss.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/06/13
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Truss 2D element********************
    template<class T>
    void Truss2D(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _E, T _A) {
        assert(_doulist.size() == 2 && _element.size() == 2);
        
        _Ke = Matrix<T>(4, 4);
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(2, std::vector<std::pair<int, int> >(2));
		_nodetoelement[0][0] = std::make_pair(_doulist[0], 0); _nodetoelement[0][1] = std::make_pair(_doulist[1], 1);
		_nodetoelement[1][0] = std::make_pair(_doulist[0], 2); _nodetoelement[1][1] = std::make_pair(_doulist[1], 3);
		
        Vector<T> dx = _x[_element[1]] - _x[_element[0]];
        T l = dx.Norm();
        T cos = dx(0)/l;
        T sin = dx(1)/l;

        T k = _E*_A/l;

        _Ke(0, 0) = k*cos*cos;  _Ke(0, 1) = k*cos*sin;  _Ke(0, 2) = -k*cos*cos; _Ke(0, 3) = -k*cos*sin;
        _Ke(1, 0) = k*cos*sin;  _Ke(1, 1) = k*sin*sin;  _Ke(1, 2) = -k*cos*sin; _Ke(1, 3) = -k*sin*sin;
        _Ke(2, 0) = -k*cos*cos; _Ke(2, 1) = -k*cos*sin; _Ke(2, 2) = k*cos*cos;  _Ke(2, 3) = k*cos*sin;
        _Ke(3, 0) = -k*cos*sin; _Ke(3, 1) = -k*sin*sin; _Ke(3, 2) = k*cos*sin;  _Ke(3, 3) = k*sin*sin;
    }
}