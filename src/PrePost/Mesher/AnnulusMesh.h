//*****************************************************************************
//  Title       :   src/PrePost/Mesher/AnnulusMesh.h
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
    //********************AnnulusMesh*******************
    template<class T>
    class AnnulusMesh{
public:
        AnnulusMesh(T _r0, T _r1, int _nr, int _nt);
        ~AnnulusMesh();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<std::vector<int> > GenerateEdges();
        template<class F>
        std::vector<std::pair<std::pair<int, int>, T> > GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond);
private:
        T r0, r1;
        int nr, nt;  
    };


    template<class T>
    AnnulusMesh<T>::AnnulusMesh(T _r0, T _r1, int _nr, int _nt){
        this->r0 = _r0;
        this->r1 = _r1;
        this->nr = _nr;
        this->nt = _nt;
    }


    template<class T>
    AnnulusMesh<T>::~AnnulusMesh(){}


    template<class T>
    std::vector<Vector<T> > AnnulusMesh<T>::GenerateNodes(){
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
    std::vector<std::vector<int> > AnnulusMesh<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nr*this->nt);
        for(int i = 0; i < this->nr; i++){
            for(int j = 0; j < this->nt; j++){
                elements[this->nt*i + j] = { this->nt*i + j%this->nt, this->nt*(i + 1) + j%this->nt, this->nt*(i + 1) + (j + 1)%this->nt, this->nt*i + (j + 1)%this->nt };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<std::vector<int> > AnnulusMesh<T>::GenerateEdges(){
        std::vector<std::vector<int> > edges = std::vector<std::vector<int> >(2*this->nt);
        for(int i = 0; i < this->nt; i++){
            edges[this->nt - i - 1] = { (i + 1)%this->nt, i };
            edges[i + this->nt] = { this->nt*this->nr + i, this->nt*this->nr + (i + 1)%this->nt };
        }
        return edges;
    }


    template<class T>
    template<class F>
    std::vector<std::pair<std::pair<int, int>, T> > AnnulusMesh<T>::GenerateFixedlist(std::vector<int> _ulist, F _iscorrespond) {
        assert(0 <= *std::min_element(_ulist.begin(), _ulist.end()));
        std::vector<std::pair<std::pair<int, int>, T> > ufixed;
        for(int i = 0; i < this->nr + 1; i++){
            T r = (this->r1 - this->r0)*i/(T)this->nr + this->r0;
            for(int j = 0; j < this->nt; j++){
                T theta = 2.0*M_PI*j/(T)this->nt;
                if(_iscorrespond(Vector<T>({ r*cos(theta), r*sin(theta) }))) {
                    for(auto ui : _ulist) {
                        ufixed.push_back({ { this->nt*i + j, ui }, T() });
                    }
                }
            }
        }
        return ufixed;
    }
}