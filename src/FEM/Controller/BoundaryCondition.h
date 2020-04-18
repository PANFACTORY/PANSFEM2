//*****************************************************************************
//  Title		:   PANSFEM2/FEM/Controller/BoundaryCondition.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/04/16
//  Copyright	:   (C)2020 TanabeYuta
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


    //********************Set Periodic boundary conditions********************
    int SetPeriodic(std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::pair<int, int> >& _ufixed) {
        for(auto nodes : _ufixed) {
            for(auto& nodetoglobali : _nodetoglobal[nodes.second]) {
                nodetoglobali = -1;
            }
        }

        int KDEGREE = 0;
        for(auto& node : _nodetoglobal) {
            for(auto& dou : node) {
                if(dou != -1) {
                    dou = KDEGREE;
                    KDEGREE++;
                }
            }
        }

        for(auto pairset : _ufixed) {
            for(int i = 0; i < _nodetoglobal[pairset.second].size(); i++) {
                _nodetoglobal[pairset.second][i] = _nodetoglobal[pairset.first][i];
            }
        }

        return KDEGREE;
    }


    //********************Remove all boundary conditions********************
    void RemoveBoundaryConditions(std::vector<std::vector<int> >& _nodetoglobal) {
        for(auto& nodetoglobali : _nodetoglobal) {
            for(auto& doui : nodetoglobali) {
                doui = 0;
            }
        }
    }
}