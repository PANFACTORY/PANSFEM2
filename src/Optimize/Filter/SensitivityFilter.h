//*****************************************************************************
//  Title       :   src/Optimize/Filter/SensitivityFilter.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/05/17
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace PANSFEM2{
    //**********Sensitivity filter class (Sigmund's scheme)**********
    template<class T>
    class SensitivityFilter{
public:
        SensitivityFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w);
        ~SensitivityFilter();
            

        std::vector<T> GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho);    //  Return filtered sensitivities

    
private:
        const int n;                                //  Number of design variables 
        std::vector<std::vector<int> > neighbors;   //  Index list of neighbor element
        std::vector<std::vector<T> > w;             //  Weight list of neighbor element    
    };


    template<class T>
    SensitivityFilter<T>::SensitivityFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w) : n(_n){
        this->neighbors = _neighbors;
        this->w = _w;
    }


    template<class T>
    SensitivityFilter<T>::~SensitivityFilter(){}


    template<class T>
    std::vector<T> SensitivityFilter<T>::GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfds){
        std::vector<T> dfds = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                dfds[i] += this->w[i][j]*_s[this->neighbors[i][j]]*_dfds[this->neighbors[i][j]];
                wsum += this->w[i][j];
            }
            dfds[i] /= wsum*_s[i];       
        }
        return dfds;
    }


    //**********Sensitivity filter class (Borrvaell's scheme)**********
    template<class T>
    class SensitivityFilter2{
public:
        SensitivityFilter2(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w);
        ~SensitivityFilter2();
            

        std::vector<T> GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho);    //  Return filtered sensitivities

    
private:
        const int n;                                //  Number of design variables 
        std::vector<std::vector<int> > neighbors;   //  Index list of neighbor element
        std::vector<std::vector<T> > w;             //  Weight list of neighbor element    
    };


    template<class T>
    SensitivityFilter2<T>::SensitivityFilter2(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w) : n(_n){
        this->neighbors = _neighbors;
        this->w = _w;
    }


    template<class T>
    SensitivityFilter2<T>::~SensitivityFilter2(){}


    template<class T>
    std::vector<T> SensitivityFilter2<T>::GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfds){
        std::vector<T> dfds = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                dfds[i] += this->w[i][j]*_s[this->neighbors[i][j]]*_dfds[this->neighbors[i][j]];
                wsum += this->w[i][j]*_s[this->neighbors[i][j]];
            }
            dfds[i] /= wsum;       
        }
        return dfds;
    }
}