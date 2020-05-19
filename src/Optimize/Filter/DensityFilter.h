//*****************************************************************************
//  Title       :   src/Optimize/Filter/DensityFilter.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/05/19
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace PANSFEM2{
    //**********Density filter class**********
    template<class T>
    class DensityFilter{
public:
        DensityFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w);
        ~DensityFilter();
            

        std::vector<T> GetFilteredVariables(std::vector<T> _s);                                 //  Return filtered variables
        std::vector<T> GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho);      //  Return filtered sensitivities

    
private:
        const int n;                                //  Number of design variables 
        std::vector<std::vector<int> > neighbors;   //  Index list of neighbor element
        std::vector<std::vector<T> > w;             //  Weight list of neighbor element    
    };


    template<class T>
    DensityFilter<T>::DensityFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w) : n(_n){
        this->neighbors = _neighbors;
        this->w = _w;
    }
    

    template<class T>
    DensityFilter<T>::~DensityFilter(){}


    template<class T>
    std::vector<T> DensityFilter<T>::GetFilteredVariables(std::vector<T> _s){
        std::vector<T> rho = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                rho[i] += this->w[i][j]*_s[this->neighbors[i][j]];
                wsum += this->w[i][j];
            }
            rho[i] /= wsum;
        }
        return rho;
    }


    template<class T>
    std::vector<T> DensityFilter<T>::GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho){
        std::vector<T> dfds = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                dfds[i] += _dfdrho[this->neighbors[i][j]]*this->w[i][j];
                wsum += this->w[i][j];
            }
            dfds[i] /= wsum;       
        }
        return dfds;
    }
}