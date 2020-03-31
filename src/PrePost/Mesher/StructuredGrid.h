//*****************************************************************************
//  Title       :   src/PrePost/Mesher/StructuredGrid.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/03/31
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#include <vector>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************Structured Grid*******************
    template<class T>
    class StructuredGrid{
public:
        StructuredGrid(T _x, T _y, int _nx, int _ny);
        ~StructuredGrid();


        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
        std::vector<int> GenerateFields(int _nu);
private:
        T x, y;
        int nx, ny;  
    };


    template<class T>
    StructuredGrid<T>::StructuredGrid(T _x, T _y, int _nx, int _ny){
        this->x = _x;
        this->y = _y;
        this->nx = _nx;
        this->ny = _ny;
    }


    template<class T>
    StructuredGrid<T>::~StructuredGrid(){}


    template<class T>
    std::vector<Vector<T> > StructuredGrid<T>::GenerateNodes(){
        std::vector<Vector<T> > nodes = std::vector<Vector<T> >((this->nx + 1)*(this->ny + 1));
        for(int i = 0; i < this->nx + 1; i++){
            for(int j = 0; j < this->ny + 1; j++){
                nodes[(this->ny + 1)*i + j] = { this->x*(i/(T)this->nx), this->y*(j/(T)this->ny) };
            }
        }
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > StructuredGrid<T>::GenerateElements(){
        std::vector<std::vector<int> > elements = std::vector<std::vector<int> >(this->nx*this->ny);
        for(int i = 0; i < this->nx; i++){
            for(int j = 0; j < this->ny; j++){
                elements[this->ny*i + j] = { (this->ny + 1)*i + j, (this->ny + 1)*(i + 1) + j, (this->ny + 1)*(i + 1) + (j + 1), (this->ny + 1)*i + (j + 1) };
            }
        }
        return elements;
    }


    template<class T>
    std::vector<int> StructuredGrid<T>::GenerateFields(int _nu){
        std::vector<int> fields = std::vector<int>((this->nx + 1)*(this->ny + 1) + 1, 0);
        for(int i = 1; i < (this->nx + 1)*(this->ny + 1) + 1; i++){
            fields[i] = fields[i - 1] + _nu;
        }
        return fields;
    }
}