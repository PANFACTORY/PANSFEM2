//*****************************************************************************
//  Title       :   src/PrePost/Mesher/SquareMesh.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/03/31
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#include <vector>
#include <algorithm>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************SquareMesh*******************
    template<class T>
    class SquareMesh{
public:
        SquareMesh(T _x, T _y, int _nx, int _ny);
        ~SquareMesh();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<std::vector<int> > GenerateEdges();
        std::vector<int> GenerateFields(int _nu);
        template<class F>
        std::vector<int> GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond);
private:
        T x, y;
        int nx, ny;  
    };


    template<class T>
    SquareMesh<T>::SquareMesh(T _x, T _y, int _nx, int _ny){
        this->x = _x;
        this->y = _y;
        this->nx = _nx;
        this->ny = _ny;
    }


    template<class T>
    SquareMesh<T>::~SquareMesh(){}


    template<class T>
    std::vector<Vector<T> > SquareMesh<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >((this->nx + 1)*(this->ny + 1));
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                nodes[(this->ny + 1)*i + j] = { this->x*(i/(T)this->nx), this->y*(j/(T)this->ny) };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > SquareMesh<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nx*this->ny);
        for(int i = 0; i < this->nx; i++){
            for(int j = 0; j < this->ny; j++){
                elements[this->ny*i + j] = { (this->ny + 1)*i + j, (this->ny + 1)*(i + 1) + j, (this->ny + 1)*(i + 1) + (j + 1), (this->ny + 1)*i + (j + 1) };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<std::vector<int> > SquareMesh<T>::GenerateEdges() {
        std::vector<std::vector<int> > edges = std::vector<std::vector<int> >(2*(this->nx + this->ny));
        for(int i = 0; i < this->nx; i++) {
            edges[i] = { (this->ny + 1)*i, (this->ny + 1)*(i + 1) };
            edges[2*this->nx + this->ny - i - 1] = { (this->ny + 1)*(i + 1) + this->ny, (this->ny + 1)*i  + this->ny };
        }
        for(int j = 0; j < this->ny; j++) {
            edges[j + this->nx] = { (this->ny + 1)*this->nx + j, (this->ny + 1)*this->nx + j + 1 };
            edges[2*(this->nx + this->ny) - j - 1] = { j + 1, j };
        }
        return edges;
    }


    template<class T>
    std::vector<int> SquareMesh<T>::GenerateFields(int _nu){
        std::vector<int> fields = std::vector<int>((this->nx + 1)*(this->ny + 1) + 1, 0);
        for(int i = 1; i < (this->nx + 1)*(this->ny + 1) + 1; i++){
            fields[i] = fields[i - 1] + _nu;
        }
        return fields;
    }


    template<class T>
    template<class F>
    std::vector<int> SquareMesh<T>::GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()) && *std::max_element(_ulist.begin(), _ulist.end()) < _nu);
        std::vector<int> isfixed;
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                if(_iscorrespond(Vector<T>({ this->x*(i/(T)this->nx), this->y*(j/(T)this->ny) }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*((this->ny + 1)*i + j) + ui);
                    }
                }
            }
        }
        return isfixed;
    }
}