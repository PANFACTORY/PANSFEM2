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
        void UpdateVariables(std::vector<T>& _xk, T _objective, std::vector<T> _dobjective, Vector<T> _constraints, std::vector<Vector<T> > _dconstraints);
        

private:
        const int n;            //Number of design variables
        const int m;            //Number of constraints


        int k;                  //Count of iteration  


        T epsf0;				//Self epsilon for objective function value
        T beforef0;             //Objective function value of before step


        T s;					//
        std::vector<T> L;       //Parameter for MMA
        std::vector<T> U;       //Parameter for MMA


        std::vector<T> xkm1;    //Design variables of previous step
        std::vector<T> xkm2;    //Design variables of previous 2 step


		T h;					//Constant for line search
		T epsalpha;				//Self epsilon for line search


		T Wy(T _r0, const Vector<T>& _rs, 
			const std::vector<T>& _p0, const std::vector<Vector<T> >& _ps, 
			const std::vector<T>& _q0, const std::vector<Vector<T> >& _qs, 
			const Vector<T>& _y, const std::vector<T>& _x);						//Function value of W(y)
		Vector<T> dWy(const Vector<T>& _rs, const std::vector<Vector<T> >& _ps,  const std::vector<Vector<T> >& _qs, const std::vector<T>& _x);		//Derivatives of W(y) 
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


		this->h = 0.01;
		
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
    void MMA<T>::UpdateVariables(std::vector<T>& _xk, T _objective, std::vector<T> _dobjective, Vector<T> _constraints, std::vector<Vector<T> > _dconstraints){
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
        T r0 = _objective;                                                                      //Objective function value at xk
        std::vector<T> p0 = std::vector<T>(this->n, T());                                       //Positive sensitivities of objective function
        std::vector<T> q0 = std::vector<T>(this->n, T());                                       //Negative sensitivities of objective function
		Vector<T> rs = _constraints;                                                 		    //Constraint function values at xk
		std::vector<Vector<T> > ps = std::vector<Vector<T> >(this->n, Vector<T>(this->m));      //Positive sensitivities of objective function
        std::vector<Vector<T> > qs = std::vector<Vector<T> >(this->n, Vector<T>(this->m));      //Negative sensitivities of objective function      
		for(int j = 0; j < this->n; j++){
            //.....Objective function.....
			if (_dobjective[j] > T()) {
				p0[j] = pow(this->U[j] - _xk[j], 2.0)*_dobjective[j];
				r0 -= p0[j] / (this->U[j] - _xk[j]);
			} else {
				q0[j] = -pow(_xk[j] - this->L[j], 2.0)*_dobjective[j];
				r0 -= q0[j] / (_xk[j] - this->L[j]);
			}

            //.....Constraint functions.....
            for(int i = 0; i < this->m; i++){
                if (_dconstraints[j](i) > T()) {
					ps[j](i) = pow(this->U[j] - _xk[j], 2.0)*_dconstraints[j](i);
					rs(i) -= ps[j](i) / (this->U[j] - _xk[j]);
				} else {
					qs[j](i) = -pow(_xk[j] - this->L[j], 2.0)*_dconstraints[j](i);
					rs(i) -= qs[j](i) / (_xk[j] - this->L[j]);
				}
            }
        }

        //----------Set movelimit----------
        std::vector<T> xmin = std::vector<T>(this->n);
		std::vector<T> xmax = std::vector<T>(this->n);
		for(int j = 0; j < this->n; j++){
			xmin[j] = std::max(0.9*this->L[j] + 0.1*_xk[j], 0.0);
			xmax[j] = std::min(0.9*this->U[j] + 0.1*_xk[j], 1.0);
		}

        //----------Loop for solving subproblem----------
		std::vector<T> xkp1 = std::vector<T>(this->n);
		Vector<T> yk = Vector<T>(this->m);

		for(int j = 0; j < this->n; j++){
			T dlmin = (p0[j] + yk*ps[j]) / pow(this->U[j] - xmin[j], 2.0) - (q0[j] + yk*qs[j]) / pow(xmin[j] - this->L[j], 2.0);
			T dlmax = (p0[j] + yk*ps[j]) / pow(this->U[j] - xmax[j], 2.0) - (q0[j] + yk*qs[j]) / pow(xmax[j] - this->L[j], 2.0);

			if (dlmin >= T()) {
				xkp1[j] = xmin[j];
			} else if (dlmax <= T()) {
				xkp1[j] = xmax[j];
			} else {
				xkp1[j] = (sqrt(p0[j] + yk*ps[j])*this->L[j] + sqrt(q0[j] + yk*qs[j])*this->U[j]) / (sqrt(p0[j] + yk*ps[j]) + sqrt(q0[j] + yk*qs[j]));
			}
		}

		Vector<T> rk = dWy(rs, ps, qs, xkp1);
		Vector<T> pk = rk;

		for(int t = 0; t < 100; t++){
			//.....Get x(y).....
			for(int j = 0; j < this->n; j++){
				T dlmin = (p0[j] + yk*ps[j]) / pow(this->U[j] - xmin[j], 2.0) - (q0[j] + yk*qs[j]) / pow(xmin[j] - this->L[j], 2.0);
				T dlmax = (p0[j] + yk*ps[j]) / pow(this->U[j] - xmax[j], 2.0) - (q0[j] + yk*qs[j]) / pow(xmax[j] - this->L[j], 2.0);

				if (dlmin >= T()) {
					xkp1[j] = xmin[j];
				} else if (dlmax <= T()) {
					xkp1[j] = xmax[j];
				} else {
					xkp1[j] = (sqrt(p0[j] + yk*ps[j])*this->L[j] + sqrt(q0[j] + yk*qs[j])*this->U[j]) / (sqrt(p0[j] + yk*ps[j]) + sqrt(q0[j] + yk*qs[j]));
				}
			}

			//.....Line search i.....
			T alpha = T(), falpha = this->Wy(r0, rs, p0, ps, q0, qs, yk + alpha*pk, xkp1);
			for(int l = 0; l < 100; l++){
				alpha += this->h;
				T falphap1 = this->Wy(r0, rs, p0, ps, q0, qs, yk + alpha*pk, xkp1);
				if(falphap1 >= falpha){
					break;
				}
				falpha = falphap1;
			}

			//.....Line serch ii.....
			T alpha1, alpha2;
			while((alpha2 - alpha1) / (alpha2 + alpha1) > epsalpha){
				
			}

			//.....Update y.....
			Vector<T> ykp1 = yk + alpha*pk;
			Vector<T> rkp1 = dWy(rs, ps, qs, xkp1);
			T beta = (rkp1*rkp1) / (rk*rk);
			Vector<T> pkp1 = rkp1 + beta*pk;

			yk = ykp1;
			rk = rkp1;
			pk = pkp1;

			if(rk.Norm() < 1.0e-5){
				break;
			}

			std::cout << yk;
		}

        //----------Update step----------
        this->k++;
    }


	template<class T>
	inline T MMA<T>::Wy(T _r0, const Vector<T>& _rs, 
			const std::vector<T>& _p0, const std::vector<Vector<T> >& _ps, 
			const std::vector<T>& _q0, const std::vector<Vector<T> >& _qs, 
			const Vector<T>& _y, const std::vector<T>& _x) {
		T value = _r0 + _y*_rs;
		for(int j = 0; j < this->n; j++){
			value += (_p0[j] + _y*_ps[j]) / (this->U[j] - _x[j]) + (_q0[j] + _y*_qs[j]) / (_x[j] - this->L[j]);
		}
		return value;
	}


	template<class T>
	inline Vector<T> dWy(const Vector<T>& _rs, const std::vector<Vector<T> >& _ps,  const std::vector<Vector<T> >& _qs, const std::vector<T>& _x) {
		Vector<T> vec = _rs;
		for(int j = 0; j < this->n; j++){
			vec += _ps[j] / (this->U[j] - _x[j]) + _qs[j] / (_x[j] - this->L[j]);
		}
		return vec;
	}
}