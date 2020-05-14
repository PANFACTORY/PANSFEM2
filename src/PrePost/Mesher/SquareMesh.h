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
        template<class F>
        std::vector<std::pair<std::pair<int, int>, T> > GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond);
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
    template<class F>
    std::vector<std::pair<std::pair<int, int>, T> > SquareMesh<T>::GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond) {
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


    //********************SquareMesh2*******************
    template<class T>
    class SquareMesh2{
public:
        SquareMesh2(T _x, T _y, int _nx, int _ny, T _rx, T _ry);
        ~SquareMesh2();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<std::vector<int> > GenerateEdges();
        template<class F>
        std::vector<std::pair<std::pair<int, int>, T> > GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond);
private:
        T x, y, rx, ry;
        int nx, ny;  
    };


    template<class T>
    SquareMesh2<T>::SquareMesh2(T _x, T _y, int _nx, int _ny, T _rx, T _ry){
        this->x = _x;
        this->y = _y;
        this->nx = _nx;
        this->ny = _ny;
        this->rx = _rx;
        this->ry = _ry;
    }


    template<class T>
    SquareMesh2<T>::~SquareMesh2(){}


    template<class T>
    std::vector<Vector<T> > SquareMesh2<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >((this->nx + 1)*(this->ny + 1));
        T rx0 = 0.5*this->x/(pow(this->rx, this->nx/2.0) - 1.0);
        T ry0 = 0.5*this->y/(pow(this->ry, this->ny/2.0) - 1.0);
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                if(i < this->nx/2.0) {
                    if(j < this->ny/2.0) {
                        nodes[(this->ny + 1)*i + j] = { rx0*(pow(this->rx, i) - 1.0), ry0*(pow(this->ry, j) - 1.0) };
                    } else if(this->ny/2.0 < j){
                        nodes[(this->ny + 1)*i + j] = { rx0*(pow(this->rx, i) - 1.0), this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) };
                    } else {
                        nodes[(this->ny + 1)*i + j] = { rx0*(pow(this->rx, i) - 1.0), 0.5*this->y };
                    }
                } else if(this->nx/2.0 < i){
                    if(j < this->ny/2.0) {
                        nodes[(this->ny + 1)*i + j] = { this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), ry0*(pow(this->ry, j) - 1.0) };
                    } else if(this->ny/2.0 < j){
                        nodes[(this->ny + 1)*i + j] = { this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) };
                    } else {
                        nodes[(this->ny + 1)*i + j] = { this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), 0.5*this->y };
                    }
                } else {
                    if(j < this->ny/2.0) {
                        nodes[(this->ny + 1)*i + j] = { 0.5*this->x, ry0*(pow(this->ry, j) - 1.0) };
                    } else if(this->ny/2.0 < j){
                        nodes[(this->ny + 1)*i + j] = { 0.5*this->x, this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) };
                    } else {
                        nodes[(this->ny + 1)*i + j] = { 0.5*this->x, 0.5*this->y };
                    }
                }
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > SquareMesh2<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nx*this->ny);
        for(int i = 0; i < this->nx; i++){
            for(int j = 0; j < this->ny; j++){
                elements[this->ny*i + j] = { (this->ny + 1)*i + j, (this->ny + 1)*(i + 1) + j, (this->ny + 1)*(i + 1) + (j + 1), (this->ny + 1)*i + (j + 1) };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<std::vector<int> > SquareMesh2<T>::GenerateEdges() {
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
    template<class F>
    std::vector<std::pair<std::pair<int, int>, T> > SquareMesh2<T>::GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()));
        std::vector<std::pair<std::pair<int, int>, T> > ufixed;
        T rx0 = 0.5*this->x/(pow(this->rx, this->nx/2.0) - 1.0);
        T ry0 = 0.5*this->y/(pow(this->ry, this->ny/2.0) - 1.0);
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                if(i < this->nx/2.0) {
                    if(j < this->ny/2.0) {
                        if(_iscorrespond(Vector<T>({ rx0*(pow(this->rx, i) - 1.0), ry0*(pow(this->ry, j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else if(this->ny/2.0 < j){
                        if(_iscorrespond(Vector<T>({ rx0*(pow(this->rx, i) - 1.0), this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else {
                        if(_iscorrespond(Vector<T>({ rx0*(pow(this->rx, i) - 1.0), 0.5*this->y }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    }
                } else if(this->nx/2.0 < i){
                    if(j < this->ny/2.0) {
                        if(_iscorrespond(Vector<T>({ this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), ry0*(pow(this->ry, j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else if(this->ny/2.0 < j){
                        if(_iscorrespond(Vector<T>({ this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else {
                        if(_iscorrespond(Vector<T>({ this->x - rx0*(pow(this->rx, this->nx - i) - 1.0), 0.5*this->y }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    }
                } else {
                    if(j < this->ny/2.0) {
                        if(_iscorrespond(Vector<T>({ 0.5*this->x, ry0*(pow(this->ry, j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else if(this->ny/2.0 < j){
                        if(_iscorrespond(Vector<T>({ 0.5*this->x, this->y - ry0*(pow(this->ry, this->ny - j) - 1.0) }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    } else {
                        if(_iscorrespond(Vector<T>({ 0.5*this->x, 0.5*this->y }))) {
                            for(auto ui : _ulist) {
                                ufixed.push_back({ { (this->ny + 1)*i + j, ui }, T() });
                            }
                        }
                    }
                }
            }
        }
        return ufixed;
    }
}