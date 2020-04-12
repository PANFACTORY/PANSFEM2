//*****************************************************************************
//  Title       :   src/PrePost/Mesher/SquareCircleAnnulus.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/04/11
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************SquareCircleAnnulus mesh*******************
    template<class T>
    class SquareCircleAnnulus{
public:
        SquareCircleAnnulus(T _a, T _b, T _r, T _p, int _nx, int _ny, int _nr);
        ~SquareCircleAnnulus();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<std::vector<int> > GenerateEdges();
        std::vector<int> GenerateFields(int _nu);
        template<class F>
        std::vector<int> GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond);
private:
        T a, b, r, p;
        int nx, ny, nr, nxy;  
    };


    template<class T>
    SquareCircleAnnulus<T>::SquareCircleAnnulus(T _a, T _b, T _r, T _p, int _nx, int _ny, int _nr){
        this->a = _a;
        this->b = _b;
        this->r = _r;
        this->p = _p;
        this->nx = _nx;
        this->ny = _ny;
        this->nr = _nr;
        this->nxy = 2*(this->nx + this->ny);
    }


    template<class T>
    SquareCircleAnnulus<T>::~SquareCircleAnnulus(){}


    template<class T>
    std::vector<Vector<T> > SquareCircleAnnulus<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >(this->nxy*(this->nr + 1));
        for(int i = 0; i < this->nr + 1; i++){
            T t = pow(i/(T)this->nr, this->p);
            for(int j = 0; j < this->ny; j++){
                T theta = 2*M_PI*(j - 0.5*this->ny)/(T)this->nxy;
                nodes[this->nxy*i + j] = { (1 - t)*this->r*cos(theta) + t*0.5*this->a, (1 - t)*this->r*sin(theta) + t*this->b*(j/(T)this->ny - 0.5) };
            }
            for(int j = 0; j < this->nx; j++){
                T theta = 2*M_PI*(j + 0.5*this->ny)/(T)this->nxy;
                nodes[this->nxy*i + j + this->ny] = { (1 - t)*this->r*cos(theta) + t*this->a*(0.5 - j/(T)this->nx), (1 - t)*this->r*sin(theta) + t*0.5*this->b };
            }
            for(int j = 0; j < this->ny; j++){
                T theta = 2*M_PI*(j + 0.5*this->ny + this->nx)/(T)this->nxy;
                nodes[this->nxy*i + j + this->ny + this->nx] = { (1 - t)*this->r*cos(theta) - t*0.5*this->a, (1 - t)*this->r*sin(theta) + t*this->b*(0.5 - j/(T)this->ny) };
            }
            for(int j = 0; j < this->nx; j++){
                T theta = 2*M_PI*(j + 1.5*this->ny + this->nx)/(T)this->nxy;
                nodes[this->nxy*i + j + 2*this->ny + this->nx] = { (1 - t)*this->r*cos(theta) + t*this->a*(j/(T)this->nx - 0.5), (1 - t)*this->r*sin(theta) - t*0.5*this->b };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > SquareCircleAnnulus<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nxy*this->nr);
        for(int i = 0; i < this->nr; i++){
            for(int j = 0; j < this->nxy; j++){
                elements[2*(this->nx + this->ny)*i + j] = { this->nxy*i + j%this->nxy, this->nxy*(i + 1) + j%this->nxy, this->nxy*(i + 1) + (j + 1)%this->nxy, this->nxy*i + (j + 1)%this->nxy };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<std::vector<int> > SquareCircleAnnulus<T>::GenerateEdges(){
        std::vector<std::vector<int> > edges = std::vector<std::vector<int> >(2*this->nxy);
        for(int i = 0; i < this->nxy; i++){
            edges[this->nxy - i - 1] = { (i + 1)%this->nxy, i };
            edges[i + this->nxy] = { this->nxy*this->nr + i, this->nxy*this->nr + (i + 1)%this->nxy };
        }
        return edges;
    }


    template<class T>
    std::vector<int> SquareCircleAnnulus<T>::GenerateFields(int _nu){
        std::vector<int> fields = std::vector<int>(this->nxy*(this->nr + 1) + 1, 0);
        for(int i = 1; i < this->nxy*(this->nr + 1) + 1; i++){
            fields[i] = fields[i - 1] + _nu;
        }
        return fields;
    }


    template<class T>
    template<class F>
    std::vector<int> SquareCircleAnnulus<T>::GenerateFixedlist(int _nu, std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()) && *std::max_element(_ulist.begin(), _ulist.end()) < _nu);
        std::vector<int> isfixed;
        for(int i = 0; i < this->nr + 1; i++){
            T t = pow(i/(T)this->nr, this->p);
            for(int j = 0; j < this->ny; j++){
                T theta = 2*M_PI*(j - 0.5*this->ny)/(T)this->nxy;
                if(_iscorrespond(Vector<T>({ (1 - t)*this->r*cos(theta) + t*0.5*this->a, (1 - t)*this->r*sin(theta) + t*this->b*(j/(T)this->ny - 0.5) }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*(this->nxy*i + j) + ui);
                    }
                }
            }
            for(int j = 0; j < this->nx; j++){
                T theta = 2*M_PI*(j + 0.5*this->ny)/(T)this->nxy;
                if(_iscorrespond(Vector<T>({ (1 - t)*this->r*cos(theta) + t*this->a*(0.5 - j/(T)this->nx), (1 - t)*this->r*sin(theta) + t*0.5*this->b }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*(this->nxy*i + j + this->ny) + ui);
                    }
                }
            }
            for(int j = 0; j < this->ny; j++){
                T theta = 2*M_PI*(j + 0.5*this->ny + this->nx)/(T)this->nxy;
                if(_iscorrespond(Vector<T>({ (1 - t)*this->r*cos(theta) - t*0.5*this->a, (1 - t)*this->r*sin(theta) + t*this->b*(0.5 - j/(T)this->ny) }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*(this->nxy*i + j + this->ny + this->nx) + ui);
                    }
                }
            }
            for(int j = 0; j < this->nx; j++){
                T theta = 2*M_PI*(j + 1.5*this->ny + this->nx)/(T)this->nxy;
                if(_iscorrespond(Vector<T>({ (1 - t)*this->r*cos(theta) + t*this->a*(j/(T)this->nx - 0.5), (1 - t)*this->r*sin(theta) - t*0.5*this->b }))) {
                    for(auto ui : _ulist) {
                        isfixed.push_back(_nu*(this->nxy*i + j + 2*this->ny + this->nx) + ui);
                    }
                }
            }
        }

        return isfixed;
    }
}