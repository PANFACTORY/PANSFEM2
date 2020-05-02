//*****************************************************************************
//  Title       :   src/PrePost/Mesher/SquareAnnulusMesh.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/04/10
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************SquareAnnulusMesh*******************
    template<class T>
    class SquareAnnulusMesh{
public:
        SquareAnnulusMesh(T _a, T _b, T _c, T _d, int _nx, int _ny, int _nt);
        ~SquareAnnulusMesh();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<std::vector<int> > GenerateEdges();
        template<class F>
        std::vector<std::pair<std::pair<int, int>, T> > GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond);
private:
        T a, b, c, d;
        int nx, ny, nt, nxy;  
    };


    template<class T>
    SquareAnnulusMesh<T>::SquareAnnulusMesh(T _a, T _b, T _c, T _d, int _nx, int _ny, int _nt){
        this->a = _a;
        this->b = _b;
        this->c = _c;
        this->d = _d;
        this->nx = _nx;
        this->ny = _ny;
        this->nt = _nt;
        this->nxy = 2*(this->nx + this->ny);
    }


    template<class T>
    SquareAnnulusMesh<T>::~SquareAnnulusMesh(){}


    template<class T>
    std::vector<Vector<T> > SquareAnnulusMesh<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >(this->nxy*(this->nt + 1));
        for(int i = 0; i < this->nt + 1; i++){
            T t = i/(T)this->nt;
            for(int j = 0; j < this->ny; j++){
                nodes[this->nxy*i + j] = { t*0.5*this->a + (1 - t)*0.5*this->c, t*this->b*(j/(T)this->ny - 0.5) + (1 - t)*this->d*(j/(T)this->ny - 0.5) };
                nodes[this->nxy*i + j + this->ny + this->nx] = { -t*0.5*this->a - (1 - t)*0.5*this->c, t*this->b*(0.5 - j/(T)this->ny) + (1 - t)*this->d*(0.5 - j/(T)this->ny) };
            }
            for(int j = 0; j < this->nx; j++){
                nodes[this->nxy*i + j + this->ny] = { t*this->a*(0.5 - j/(T)this->nx) + (1 - t)*this->c*(0.5 - j/(T)this->nx), t*0.5*this->b + (1 - t)*0.5*this->d };
                nodes[this->nxy*i + j + 2*this->ny + this->nx] = { t*this->a*(j/(T)this->nx - 0.5) + (1 - t)*this->c*(j/(T)this->nx - 0.5), -t*0.5*this->b - (1 - t)*0.5*this->d };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > SquareAnnulusMesh<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nxy*this->nt);
        for(int i = 0; i < this->nt; i++){
            for(int j = 0; j < this->nxy; j++){
                elements[2*(this->nx + this->ny)*i + j] = { this->nxy*i + j%this->nxy, this->nxy*(i + 1) + j%this->nxy, this->nxy*(i + 1) + (j + 1)%this->nxy, this->nxy*i + (j + 1)%this->nxy };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<std::vector<int> > SquareAnnulusMesh<T>::GenerateEdges(){
        std::vector<std::vector<int> > edges = std::vector<std::vector<int> >(2*this->nxy);
        for(int i = 0; i < this->nxy; i++){
            edges[this->nxy - i - 1] = { (i + 1)%this->nxy, i };
            edges[i + this->nxy] = { this->nxy*this->nt + i, this->nxy*this->nt + (i + 1)%this->nxy };
        }
        return edges;
    }


    template<class T>
    template<class F>
    std::vector<std::pair<std::pair<int, int>, T> > SquareAnnulusMesh<T>::GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()));
        std::vector<std::pair<std::pair<int, int>, T> > ufixed;
        for(int i = 0; i < this->nt + 1; i++){
            T t = i/(T)this->nt;
            for(int j = 0; j < this->ny; j++){
                if(_iscorrespond(Vector<T>({ (1 - t)*0.5*this->a + t*0.5*this->c, (1 - t)*this->b*(j/(T)this->ny - 0.5) + t*this->d*(j/(T)this->ny - 0.5) }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { this->nxy*i + j, ui }, T() });
                    }
                }
            }
            for(int j = 0; j < this->nx; j++){
                if(_iscorrespond(Vector<T>({ (1 - t)*this->a*(0.5 - j/(T)this->nx) + t*this->c*(0.5 - j/(T)this->nx), (1 - t)*0.5*this->b + t*0.5*this->d }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { this->nxy*i + j + this->ny, ui }, T() });
                    }
                }
            }
            for(int j = 0; j < this->ny; j++){
                if(_iscorrespond(Vector<T>({ -(1 - t)*0.5*this->a - t*0.5*this->c, (1 - t)*this->b*(0.5 - j/(T)this->ny) + t*this->d*(0.5 - j/(T)this->ny) }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { this->nxy*i + j + this->ny + this->nx, ui }, T() });
                    }
                }
            }
            for(int j = 0; j < this->nx; j++){
                if(_iscorrespond(Vector<T>({ (1 - t)*this->a*(j/(T)this->nx - 0.5) + t*this->c*(j/(T)this->nx - 0.5), -(1 - t)*0.5*this->b - t*0.5*this->d }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { this->nxy*i + j + 2*this->ny + this->nx, ui }, T() });
                    }
                }
            }
        }
        return ufixed;
    }
}