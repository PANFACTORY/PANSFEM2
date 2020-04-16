//*****************************************************************************
//Title		:src/FEM/Controller/Assembling2.h
//Author	:Tanabe Yuta
//Date		:2020/04/16
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/LILCSR.h"
#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Assembling global matrix********************
    template<class T>
    void Assembling(LILCSR<T>& _K, std::vector<T>& _F, std::vector<Vector<T> >& _u, Matrix<T>& _Ke, const std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element) {
        for(int i = 0; i < _element.size(); i++) {
            for(auto doui : _nodetoelement[i]) {
                if(_nodetoglobal[_element[i]][doui.first] != -1) {
                    for(int j = 0; j < _element.size(); j++) {
                        for(auto douj : _nodetoelement[j]) {
                            //----------Dirichlet condition NOT imposed----------
                            if(_nodetoglobal[_element[j]][douj.first] != -1) {
                                _K.set(_nodetoglobal[_element[i]][doui.first], _nodetoglobal[_element[j]][douj.first], _K.get(_nodetoglobal[_element[i]][doui.first], _nodetoglobal[_element[j]][douj.first]) + _Ke(doui.second, douj.second));
                            }
                            //----------Dirichlet condition imposed----------
                            else {
                                _F[_nodetoglobal[_element[i]][doui.first]] -= _Ke(doui.second, douj.second)*_u[_element[j]](douj.first);
                            }
                        }
                    }
                }
            }
        }
    }


    //********************Assembling global matrix********************
    template<class T>
    void Assembling(LILCSR<T>& _K, Matrix<T>& _Ke, const std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element) {
        for(int i = 0; i < _element.size(); i++) {
            for(auto doui : _nodetoelement[i]) {
                if(_nodetoglobal[_element[i]][doui.first] != -1) {
                    for(int j = 0; j < _element.size(); j++) {
                        for(auto douj : _nodetoelement[j]) {
                            //----------Dirichlet condition NOT imposed----------
                            if(_nodetoglobal[_element[j]][douj.first] != -1) {
                                _K.set(_nodetoglobal[_element[i]][doui.first], _nodetoglobal[_element[j]][douj.first], _K.get(_nodetoglobal[_element[i]][doui.first], _nodetoglobal[_element[j]][douj.first]) + _Ke(doui.second, douj.second));
                            }
                        }
                    }
                }
            }
        }
    }


    //********************Assembling global vector********************
    template<class T>
    void Assembling(std::vector<T>& _F, Vector<T>& _Fe, const std::vector<std::vector<int> >& _nodetoglobal, const std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element) {
        for(int i = 0; i < _element.size(); i++) {
            for(auto doui : _nodetoelement[i]) {
                if(_nodetoglobal[_element[i]][doui.first] != -1) {
                    _F[_nodetoglobal[_element[i]][doui.first]] += _Fe(doui.second);
                }
            }
        }
    }


    //********************Assembling Neumann boundary conditions********************
    template<class T>
    void Assembling(std::vector<T>& _F, const std::vector<std::pair<std::pair<int, int>, T> >& _f, const std::vector<std::vector<int> >& _nodetoglobal) {
        for(auto doui : _f) {
            if(_nodetoglobal[doui.first.first][doui.first.second] != -1) {
                _F[_nodetoglobal[doui.first.first][doui.first.second]] += doui.second;
            }
        }
    }


    //********************Disassembling global result********************
    template<class T>
    void Disassembling(std::vector<Vector<T> >& _u, const std::vector<T>& _result, const std::vector<std::vector<int> >& _nodetoglobal) {
        for(int i = 0; i < _nodetoglobal.size(); i++) {
            for(int j = 0; j < _nodetoglobal[i].size(); j++) {
                if(_nodetoglobal[i][j] != -1) {
                    _u[i](j) += _result[_nodetoglobal[i][j]]; 
                }
            }
        }
    }


    //********************Renumbering********************
    int Renumbering(std::vector<std::vector<int> >& _nodetoglobal) {
        int KDEGREE = 0;
        for(auto& node : _nodetoglobal) {
            for(auto& dou : node) {
                if(dou != -1) {
                    dou = KDEGREE;
                    KDEGREE++;
                }
            }
        }
        return KDEGREE;
    }
}