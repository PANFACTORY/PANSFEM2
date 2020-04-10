//*****************************************************************************
//  Title       :   src/PrePost/Mesher/Annulus.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/04/09
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************Annulus mesh*******************
    template<class T>
    class Annulus{
public:
        Annulus(T _r0, T _r1, int _nr, int _nt);
        ~Annulus();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<int> GenerateFields(int _nu);
        template<class F>
        std::vector<int> GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond);
private:
        T r0, r1;
        int nr, nt;  
    };


    template<class T>
    Annulus<T>::Annulus(T _r0, T _r1, int _nr, int _nt){
        this->r0 = _r0;
        this->r1 = _r1;
        this->nr = _nr;
        this->nt = _nt;
    }


    template<class T>
    Annulus<T>::~Annulus(){}


    template<class T>
    std::vector<Vector<T> > Annulus<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >((this->nr + 1)*this->nt);
        for(int i = 0; i < this->nr + 1; i++){
            T r = (this->r1 - this->r0)*i/(T)this->nr + this->r0;
            for(int j = 0; j < this->nt; j++){
                T theta = 2.0*M_PI*j/(T)this->nt;
                nodes[this->nt*i + j] = { r*cos(theta), r*sin(theta) };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > Annulus<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nr*this->nt);
        for(int i = 0; i < this->nr; i++){
            for(int j = 0; j < this->nt; j++){
                elements[this->nt*i + j] = { this->nt*i + j%this->nt, this->nt*(i + 1) + j%this->nt, this->nt*(i + 1) + (j + 1)%this->nt, this->nt*i + (j + 1)%this->nt };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<int> Annulus<T>::GenerateFields(int _nu){
        std::vector<int> fields = std::vector<int>((this->nr + 1)*this->nt + 1, 0);
        for(int i = 1; i < (this->nr + 1)*this->nt + 1; i++){
            fields[i] = fields[i - 1] + _nu;
        }
        return fields;
    }


    template<class T>
    template<class F>
    std::vector<int> Annulus<T>::GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()) && *std::max_element(_ulist.begin(), _ulist.end()) < _nu);
        std::vector<int> isfixed;
        for(int i = 0; i < this->nr + 1; i++){
            T r = (this->r1 - this->r0)*i/(T)this->nr + this->r0;
            for(int j = 0; j < this->nt; j++){
                T theta = 2.0*M_PI*j/(T)this->nt;
                if(_iscorrespond(Vector<T>({ r*cos(theta), r*sin(theta) }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*(this->nt*i + j) + ui);
                    }
                }
            }
        }
        return isfixed;
    }
}