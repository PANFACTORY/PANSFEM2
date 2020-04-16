//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/BoundaryCondition2.h
//Author	:Tanabe Yuta
//Date		:2020/04/16
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
    //********************Set Dirichlet boundary conditions********************
    template<class T>
    void SetDirichlet(std::vector<Vector<T> >& _u, std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::pair<std::pair<int, int>, T> >& _ufixed) {
        for(auto ufixed : _ufixed) {
            _u[ufixed.first.first](ufixed.first.second) = ufixed.second;
            _nodetoglobal[ufixed.first.first][ufixed.first.second] = -1;
        }
    }


    //********************Set Dirichlet boundary conditions********************
    template<class T>
    void SetDirichlet(std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::pair<std::pair<int, int>, T> >& _ufixed) {
        for(auto ufixed : _ufixed) {
            _nodetoglobal[ufixed.first.first][ufixed.first.second] = -1;
        }
    }
}