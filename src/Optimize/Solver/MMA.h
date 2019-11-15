//*****************************************************************************
//  Title       :src/Optimize/Solver/MMA.h
//  Author      :Tanabe Yuta
//  Date        :2019/11/15
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2{
    //********************Optimizational solver with MMA********************
    template<class T>
    class MMA{
public:
        MMA(int _n, int _m);
        ~MMA();


        bool IsConvergence(T _currentf0);
        void UpdateVariables(std::vector<T>& _xk, T _objective, const std::vector<T>& _dobjective, const Vector<T>& _constraints, const std::vector<Vector<T> >& _dconstraints);
        

private:
        const int n;            //Number of design variables
        const int m;            //Number of constraints


        int k;                  //Count of iteration  


        const T epsf0;          //Self epsilon for objective function value
        T beforef0;             //Objective function value of before step


        const T s;              //
        std::vector<T> L;       //Parameter for MMA
        std::vector<T> U;       //Parameter for MMA


        std::vector<T> xkm1;    //Design variables of previous step
        std::vector<T> xkm2;    //Design variables of previous 2 step
    };


    template<class T>
    MMA<T>::MMA(int _n, int _m) : n(_n), m(_m) {
        this->k = 0;


        this->epsf0 = 1.0e-5;
        this->beforef0 = T();


        this->s = 0.7;
        this->L = std::vector<T>(this->n);
	    this->U = std::vector<T>(this->n);


        this->xkm1 = std::vector<T>(this->n, T());
        this->xkm2 = std::vector<T>(this->n, T());
    }


    template<class T>
    MMA<T>::~MMA<T>(){}


    template<class T>
    bool MMA<T>::IsConvergence(T _currentf0){
        if(fabs((_currentf0 - this->beforef0) / (_currentf0 + this->beforef0)) < this->epsf0) {
			return true;
		}
        return false;
    }


    template<class T>
    void UpdateVariables(std::vector<T>& _xk, T _objective, const std::vector<T>& _dobjective, const Vector<T>& _constraints, const std::vector<Vector<T> >& _dconstraints){
        //----------Set parameter U and L----------
		if(this->k < 2){
			for(int j = 0; j < this->n; j++){
				this->L[j] = _xk[j] - (1.0 - 0.0);
				this->U[j] = _xk[j] + (1.0 - 0.0);
			}
		} else {
			for(int j = 0; j < this->n; j++){
				if((_xk[j] - this->xkm1[j])*(this->xkm1[j] - this->xkm2[j]) < 0.0){
					this->L[j] = _xk[j] - this->s*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->s*(this->U[j] - this->xkm1[j]);
				} else {
					this->L[j] = _xk[j] - (this->xkm1[j] - this->L[j])/this->s;
					this->U[j] = _xk[j] + (this->U[j] - this->xkm1[j])/this->s;
				}
			}
		}

        //----------Similerize objective function and constraint functions----------
        T r0 = T();                                                                             //Objective function value at xk
        std::vector<T> p0 = std::vector<T>(this->n, T());                                       //Positive sensitivities of objective function
        std::vector<T> q0 = std::vector<T>(this->n, T());                                       //Negative sensitivities of objective function
        Vector<T> bs = Vector<T>(this->m);                                                      //Constraint function values at xk
        std::vector<Vector<T> > ps = std::vector<Vector<T> >(this->n, Vector<T>(this->m));      //Positive sensitivities of objective function
        std::vector<Vector<T> > qs = std::vector<Vector<T> >(this->n, Vector<T>(this->m));      //Negative sensitivities of objective function
        for(int j = 0; j < this->n; j++){
            //.....Objective function.....

            //.....Constraint functions.....
            for(int i = 0; i < this->m; i++){
                
            }
        }

        //----------Set movelimit----------
        std::vector<T> xmin = std::vector<T>(this->n);
		std::vector<T> xmax = std::vector<T>(this->n);
		for(int j = 0; j < this->n; j++){
			xmin[j] = std::max(0.9*this->L[i] + 0.1*_xk[i], 0.0);
			xmax[j] = std::min(0.9*this->U[i] + 0.1*_xk[i], 1.0);
		}

        //----------Loop for solving subproblem----------
		std::vector<T> xkp1 = std::vector<T>(this->n);
		Vector<T> yk = Vector<T>(this->m);

		Vector<T> rk = -b;
		for(int j = 0; j < this->n; j++){
			rk += ps[j] / (this->U[j] - _xk[j]) + qs[j] / (_xk[j] - this->L[j]);
		} 
		Vector<double> pk = rk;

		for(int t = 0; t < 100; t++){
			//.....Get x(y).....
			for(int i = 0; i < elements.size(); i++){
				double dlmin = (p0[i] + yk*ps[i]) / pow(U[i] - xmin[i], 2.0) - (q0[i] + yk*qs[i]) / pow(xmin[i] - L[i], 2.0);
				double dlmax = (p0[i] + yk*ps[i]) / pow(U[i] - xmax[i], 2.0) - (q0[i] + yk*qs[i]) / pow(xmax[i] - L[i], 2.0);

				if (dlmin >= 0.0) {
					xkp1[i] = xmin[i];
				} else if (dlmax <= 0.0) {
					xkp1[i] = xmax[i];
				} else {
					xkp1[i] = (sqrt(p0[i] + yk*ps[i])*L[i] + sqrt(q0[i] + yk*qs[i])*U[i]) / (sqrt(p0[i] + yk*ps[i]) + sqrt(q0[i] + yk*qs[i]));
				}
			}

			//.....Get y.....
			double alpha = 1;
			Vector<double> ykp1 = yk + alpha*pk;
			Vector<double> rkp1 = -b;
			for(int i = 0; i < elements.size(); i++){
				rkp1 += ps[i] / (U[i] - xkp1[i]) + qs[i] / (xkp1[i] - L[i]);
			}
			double beta = (rkp1*rkp1) / (rk*rk);
			Vector<double> pkp1 = rkp1 + beta*pk;
			yk = ykp1;
			rk = rkp1;
			pk = pkp1;
			if(rk.Norm() < 1.0e-5){
				break;
			}

			std::cout << rk; 
		}

        //----------Update step----------
        this->k++;
    }
}