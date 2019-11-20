//*****************************************************************************
//  Title       :src/Optimize/Solver/PrimalDual.h
//  Author      :Tanabe Yuta
//  Date        :2019/11/21
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Vector.h"
#include "../../LinearAlgebra/Models/Matrix.h"


namespace PANSFEM2{
    //**********Optimizational solver with Primal-Dual Inner Point Method**********
    template<class T>
    class PDIPM{
public:
        PDIPM();
        ~PDIPM();
private:
    };


    template<class T>
    PDIPM<T>::PDIPM(){}


    template<class T>
    PDIPM<T>::~PDIPM(){}
}