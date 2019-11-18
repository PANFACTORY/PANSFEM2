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


		T alpha0;				//Constant for line search
		T rho;					//Self epsilon for line search


		T V;					//Constant for penalty


		T Wy(T _r0, Vector<T> _rs, 
			std::vector<T> _p0, std::vector<Vector<T> > _ps, 
			std::vector<T> _q0, std::vector<Vector<T> > _qs, 
			Vector<T> _y, std::vector<T> _x);						//Function value of W(y)
		Vector<T> dWy(Vector<T> _rs, std::vector<Vector<T> > _ps,  std::vector<Vector<T> > _qs, Vector<T> _y, std::vector<T> _x);		//Derivatives of W(y) 
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


		this->alpha0 = 1.0;
		this->rho = 0.5;
		

		this->V = 1.0e-3;
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
			xmin[j] = std::max(0.9*this->L[j] + 0.1*_xk[j], 1.0e-10);
			xmax[j] = std::min(0.9*this->U[j] + 0.1*_xk[j], 1.0);
		}

        //----------Loop for solving subproblem----------
		Vector<T> yl = Vector<T>(std::vector<T>(this->m, 1.0e-20));					
		Matrix<T> Bl = Identity<T>(this->m);
		for(int l = 0; l < 500; l++){
			//.....Get x(y).....
			for(int j = 0; j < this->n; j++){
				if ((p0[j] + yl*ps[j]) / pow(this->U[j] - xmin[j], 2.0) - (q0[j] + yl*qs[j]) / pow(xmin[j] - this->L[j], 2.0) >= T()) {
					_xk[j] = xmin[j];
				} else if ((p0[j] + yl*ps[j]) / pow(this->U[j] - xmax[j], 2.0) - (q0[j] + yl*qs[j]) / pow(xmax[j] - this->L[j], 2.0) <= T()) {
					_xk[j] = xmax[j];
				} else {
					_xk[j] = (sqrt(p0[j] + yl*ps[j])*this->L[j] + sqrt(q0[j] + yl*qs[j])*this->U[j]) / (sqrt(p0[j] + yl*ps[j]) + sqrt(q0[j] + yl*qs[j]));
				}
			}
			
			//.....Get objective value and constraint at l.....
			Vector<T> df = this->dWy(rs, ps, qs, yl, _xk);
			Vector<T> g = -yl;
			Matrix<T> dg = -Identity<T>(this->m);

			//.....Solve subproblem.....
			Matrix<T> A = Bl.Hstack(dg).Vstack(dg.Transpose().Hstack(Matrix<T>(this->m, this->m)));
			Vector<T> b = -df.Vstack(-g);
			Vector<T> yz = A.Inverse()*b;
			Vector<T> dyl = yz.Segment(0, this->m);
			Vector<T> zlp1 = yz.Segment(this->m, this->m*2);

			//.....Check KKT condition.....
			

			//.....Get step size with Armijo condition.....
			T alpha = this->alpha0;
			T c = 0.5;
			for(int t = 0; t < 1000000; t++){
				//.....Get x(y).....
				Vector<T> ylp1 = yl + alpha*dyl;
				std::vector<T> xkp1 = std::vector<T>(this->n);
				for(int j = 0; j < this->n; j++){
					if ((p0[j] + ylp1*ps[j]) / pow(this->U[j] - xmin[j], 2.0) - (q0[j] + ylp1*qs[j]) / pow(xmin[j] - this->L[j], 2.0) >= T()) {
						xkp1[j] = xmin[j];
					} else if ((p0[j] + ylp1*ps[j]) / pow(this->U[j] - xmax[j], 2.0) - (q0[j] + ylp1*qs[j]) / pow(xmax[j] - this->L[j], 2.0) <= T()) {
						xkp1[j] = xmax[j];
					} else {
						xkp1[j] = (sqrt(p0[j] + ylp1*ps[j])*this->L[j] + sqrt(q0[j] + ylp1*qs[j])*this->U[j]) / (sqrt(p0[j] + ylp1*ps[j]) + sqrt(q0[j] + ylp1*qs[j]));
					}
				}

				//.....Check Armijo condition.....
				if(Wy(r0, rs, p0, ps, q0, qs, ylp1, xkp1) <= Wy(r0, rs, p0, ps, q0, qs, yl, xkp1) + c*dWy(rs, ps, qs, yl, xkp1)*alpha*dyl){
					break;
				}
				alpha *= rho;
			}
			Vector<T> ylp1 = yl + alpha*dyl;
			yl = ylp1;

			if(fabs(alpha) < 1.0e-15){
				break;
			}

			//.....Update Bl with BFGS.....
			




			std::cout << alpha << "\t" << yl(0) << "\t" << yz.Transpose();	
		}
		

        //----------Update step----------
        this->k++;
		this->xkm2 = this->xkm1;
		this->xkm1 = _xk;
    }


	template<class T>
	T MMA<T>::Wy(T _r0, Vector<T> _rs, 
			std::vector<T> _p0, std::vector<Vector<T> > _ps, 
			std::vector<T> _q0, std::vector<Vector<T> > _qs, 
			Vector<T> _y, std::vector<T> _x) {
		T value = _r0 + _y*_rs;
		for(int j = 0; j < this->n; j++){
			value += (_p0[j] + _y*_ps[j]) / (this->U[j] - _x[j]) + (_q0[j] + _y*_qs[j]) / (_x[j] - this->L[j]);
		}
		return -value;
	}


	template<class T>
	Vector<T> MMA<T>::dWy(Vector<T> _rs, std::vector<Vector<T> > _ps,  std::vector<Vector<T> > _qs, Vector<T> _y, std::vector<T> _x) {
		Vector<T> vec = _rs;
		for(int j = 0; j < this->n; j++){
			vec += _ps[j] / (this->U[j] - _x[j]) + _qs[j] / (_x[j] - this->L[j]);
		}
		return -vec;
	}
}