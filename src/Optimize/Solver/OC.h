//*****************************************************************************
//  Title       :src/Optimize/Solver/OC.h
//  Author      :Tanabe Yuta
//  Date        :2020/03/18
//  Copyright   :(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace PANSFEM2{
    //********************Optimizational solver with OC********************
    template<class T>
    class OC{
public:
        OC(int _n, T _iota, T _lambdamin, T _lambdamax, T _lambdaeps, T _movelimit, const std::vector<T>& _xmin, const std::vector<T>& _xmax);
        ~OC();


        bool IsConvergence(T _currentf0);
        template<class F>
        void UpdateVariables(std::vector<T>& _xk, T _f, std::vector<T> _dfdx, T _g, std::vector<T> _dgdx, F _gkp1);


private:
        //----------Parameters for solver----------
		int k;                  //  Counter for outer loop
        const int n;            //  Number of design variables
        T previousvalue;        //  Previous function value
        T epsvalue;             //  Convergence parameter of objective value
        std::vector<T> xmin;    //  Minimum value of design variable
        std::vector<T> xmax;    //  Maximum value of design variable

        //----------Parameters for OC----------
        T iota;
        T lambdamin;
        T lambdamax;
        T lambdaeps;
        T movelimit;
    };


    template<class T>
    OC<T>::OC(int _n, T _iota, T _lambdamin, T _lambdamax, T _lambdaeps, T _movelimit, const std::vector<T>& _xmin, const std::vector<T>& _xmax) : n(_n){
        //----------Initialize solver parameter----------
        this->k = 0;
        this->previousvalue = T();
        this->epsvalue = 1.0e-5;
        this->xmin = _xmin;
        this->xmax = _xmax;

        //----------Set default OC parameters----------
        this->iota = _iota;
        this->lambdamin = _lambdamin;
        this->lambdamax = _lambdamax;
        this->lambdaeps = _lambdaeps;
        this->movelimit = _movelimit;
    }


    template<class T>
    OC<T>::~OC(){}


    template<class T>
    bool OC<T>::IsConvergence(T _currentf0){
        if(fabs(_currentf0 - this->previousvalue)/(_currentf0 + this->previousvalue) < this->epsvalue){
            return true;
        } 
        return false;
    }


    template<class T>
    template<class F>
    void OC<T>::UpdateVariables(std::vector<T>& _xk, T _f, std::vector<T> _dfdx, T _g, std::vector<T> _dgdx, F _gkp1){
        //----------Get updated design variables with OC method----------
		T lambda0 = this->lambdamin, lambda1 = this->lambdamax, lambda;
		std::vector<T> xkp1 = std::vector<T>(this->n);
		while((lambda1 - lambda0)/(lambda1 + lambda0) > this->lambdaeps){
			lambda = 0.5*(lambda1 + lambda0);

			for (int i = 0; i < this->n; i++) {
				xkp1[i] = pow(-_dfdx[i]/(_dgdx[i]*lambda), this->iota)*_xk[i];
				if(xkp1[i] < std::max(T(), (1.0 - this->movelimit)*_xk[i])) {
					xkp1[i] = std::max(T(), (1.0 - this->movelimit)*_xk[i]);
				} else if(xkp1[i] > std::min(1.0, (1.0 + this->movelimit)*_xk[i])) {
					xkp1[i] = std::min(1.0, (1.0 + this->movelimit)*_xk[i]);
				}
			}

			if (_gkp1(xkp1) > T()) {
				lambda0 = lambda;
			} else {
				lambda1 = lambda;
			}
		}

        std::cout << lambda;

		//----------Update design variables and objective function value----------
        this->previousvalue = _f;
        this->k++;
		_xk = xkp1;
    }
}