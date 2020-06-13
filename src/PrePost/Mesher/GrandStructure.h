//*****************************************************************************
//  Title       :   src/PrePost/Mesher/GrandStructure.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/06/13
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#include <vector>
#include <algorithm>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************GrandStructure2D********************
    template<class T>
    class GrandStructure2D{
public:
        GrandStructure2D(T _x, T _y, int _nx, int _ny, T _lmax);
        ~GrandStructure2D();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        template<class F>
        std::vector<std::pair<std::pair<int, int>, T> > GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond);


private:
        T x, y, lmax;
        int nx, ny;
    };


    template<class T>
    GrandStructure2D<T>::GrandStructure2D(T _x, T _y, int _nx, int _ny, T _lmax){
        this->x = _x;
        this->y = _y;
        this->nx = _nx;
        this->ny = _ny;
        this->lmax = _lmax;
    }


    template<class T>
    GrandStructure2D<T>::~GrandStructure2D(){}


    template<class T>
    std::vector<Vector<T> > GrandStructure2D<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >((this->nx + 1)*(this->ny + 1));
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                nodes[(this->ny + 1)*i + j] = { this->x*(i/(T)this->nx), this->y*(j/(T)this->ny) };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > GrandStructure2D<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >();
        for(int i = 0; i < this->nx + 1; i++) {
            for(int j = 0; j < this->ny + 1; j++) {
                for(int k = 0; k < this->nx + 1; k++) {
                    for(int l = 0; l < this->ny + 1; l++) {
                        if((this->ny + 1)*i + j > (this->ny + 1)*k + l && pow(this->x*(i - k)/(T)this->nx, 2.0) + pow(this->y*(j - l)/(T)this->ny, 2.0) < pow(this->lmax, 2.0)) {
                            elements.push_back({ (this->ny + 1)*i + j, (this->ny + 1)*k + l });
                        }
                    }
                }
            } 
        }
        return elements;
    }


    template<class T>
    template<class F>
    std::vector<std::pair<std::pair<int, int>, T> > GrandStructure2D<T>::GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()));
        std::vector<std::pair<std::pair<int, int>, T> > ufixed;
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                if(_iscorrespond(Vector<T>({ this->x*(i/(T)this->nx), this->y*(j/(T)this->ny) }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                    }
                }
            }
        }
        return ufixed;
    }
}