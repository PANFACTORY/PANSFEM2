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


        T sm;					//
		T sn;					//
		T sp;					//
        std::vector<T> L;       //Parameter for MMA
        std::vector<T> U;       //Parameter for MMA


        std::vector<T> xkm1;    //Design variables of previous step
        std::vector<T> xkm2;    //Design variables of previous 2 step


		T alpha0;				//Constant for line search
		T rho;					//Self epsilon for line search
		T c1;


		T Mc;
		T tau;
		T mu0;
		T mumin;


		Vector<T> y;


		T Wy(T _r0, Vector<T> _rs, 
			std::vector<T> _p0, std::vector<Vector<T> > _ps, 
			std::vector<T> _q0, std::vector<Vector<T> > _qs, 
			Vector<T> _y, std::vector<T> _x);						//Function value of W(y)
		Vector<T> dWy(Vector<T> _rs, std::vector<Vector<T> > _ps,  std::vector<Vector<T> > _qs, std::vector<T> _x);		//Derivatives of W(y) 
		T LogBallier(T _mu, Vector<T> _y);
	};


    template<class T>
    MMA<T>::MMA(int _n, int _m) : n(_n), m(_m) {
        this->k = 0;


        this->epsf0 = 1.0e-5;
        this->beforef0 = T();


        this->sm = 0.7;
		this->sn = 1.0;
		this->sp = 1.2;
        this->L = std::vector<T>(this->n);
	    this->U = std::vector<T>(this->n);


        this->xkm1 = std::vector<T>(this->n, T());
        this->xkm2 = std::vector<T>(this->n, T());


		this->alpha0 = 1.0;
		this->rho = 0.5;
		this->c1 = 1.0e-4;
		

		this->Mc = 1.0;
		this->tau = 0.5;
		this->mu0 = 1.0;
		this->mumin = 1.0e-10;


		this->y = Vector<T>(std::vector<T>(this->m, 1.0e10));
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
				this->L[j] = _xk[j] - 0.5*(1.0 - 0.0);
				this->U[j] = _xk[j] + 0.5*(1.0 - 0.0);
			}
		} else {
			for(int j = 0; j < this->n; j++){
				if((_xk[j] - this->xkm1[j])*(this->xkm1[j] - this->xkm2[j]) < T()){
					this->L[j] = _xk[j] - this->sm*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->sm*(this->U[j] - this->xkm1[j]);
				} else if((_xk[j] - this->xkm1[j])*(this->xkm1[j] - this->xkm2[j]) > T()){
					this->L[j] = _xk[j] - this->sp*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->sp*(this->U[j] - this->xkm1[j]);
				} else {
					this->L[j] = _xk[j] - this->sn*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->sn*(this->U[j] - this->xkm1[j]);
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

		//----------Solve subproblem with Primal-Dual Inner Point Method----------
		T mu = this->mu0;

		//----------External loop----------
		Vector<T> z = Vector<T>(std::vector<T>(this->m, 1.0));
		while(mu > this->mumin){
			//.....Internal loop.....
			Matrix<T> Bl = Identity<T>(this->m);
			for(int l = 0; l < 100; l++){
				//...Check KKT condition...
				Matrix<T> Diagy = Diagonal<T>(y);
				Matrix<T> Diagz = Diagonal<T>(z);
				Vector<T> e = Vector<T>(std::vector<T>(this->m, 1.0));
				Vector<T> dL = this->dWy(rs, ps, qs, _xk) - z;
				Vector<T> r = dL.Vstack(Diagy*Diagz*e - mu*e);
				if(r.Norm() < this->Mc*mu){
					std::cout << "!!!"; 
					break;
				}

				//...Get search direction...
				Matrix<T> A = Bl.Hstack(-Identity<T>(this->m)).Vstack(Diagz.Hstack(Diagy));
				Vector<T> dyz = -A.Inverse()*r;
				Vector<T> dy = dyz.Segment(0, this->m);
				Vector<T> dz = dyz.Segment(this->m, 2*this->m);

				//...Get step size...
				
				//	Get maximam of alphay
				T ganma = 0.9995;
				T alphay = 1.0 / ganma;	
				for(int j = 0; j < this->m; j++){
					if(dy(j) < T() && alphay > -y(j) / dy(j)){
						alphay = -y(j) / dy(j);
					}
				}
				alphay *= ganma;

				//	Check Armijo condition
				T F0 = this->Wy(r0, rs, p0, ps, q0, qs, y, _xk) + this->LogBallier(mu, y);
				T dF0 = -dy*((Bl - Diagy.Inverse()*Diagz)*dy);
				std::vector<T> xkp1 = std::vector<T>(this->n);
				for(int t = 0; t < 100; t++){
					//.....Get updated y.....
					Vector<T> ylp1 = y + alphay*dy;

					//.....Get updated x(y).....
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
					if(this->Wy(r0, rs, p0, ps, q0, qs, ylp1, xkp1) + this->LogBallier(mu, ylp1) <= F0 + this->c1*alphay*dF0){
						break;
					}

					alphay *= this->rho;
				}

				//	Get cLki and cUki
				T ML = 100.0;
				T MU = 100.0;
				std::vector<T> cL = std::vector<T>(this->m);
				std::vector<T> cU = std::vector<T>(this->m);
				for(int j = 0; j < this->m; j++){
					cL[j] = std::min(mu/ML, (y(j) + alphay*dy(j))*z(j));
					cU[j] = std::max(mu*MU, (y(j) + alphay*dy(j))*z(j));
				}

				//	Get step size alphaz
				T alphaz = 1.0;
				for(int j = 0; j < this->m; j++){
					T tmp = std::max((cL[j]/(y(j) + alphay*dy(j)) - z(j))/dz(j), (cU[j]/(y(j) + alphay*dy(j)) - z(j))/dz(j));
					if(alphaz > tmp){
						alphaz = tmp;
					}
				}

				//...Update y and z...
				y += alphay*dy;
				z += alphaz*dz;
				for(int j = 0; j < this->n; j++){
					if ((p0[j] + y*ps[j]) / pow(this->U[j] - xmin[j], 2.0) - (q0[j] + y*qs[j]) / pow(xmin[j] - this->L[j], 2.0) >= T()) {
						_xk[j] = xmin[j];
					} else if ((p0[j] + y*ps[j]) / pow(this->U[j] - xmax[j], 2.0) - (q0[j] + y*qs[j]) / pow(xmax[j] - this->L[j], 2.0) <= T()) {
						_xk[j] = xmax[j];
					} else {
						_xk[j] = (sqrt(p0[j] + y*ps[j])*this->L[j] + sqrt(q0[j] + y*qs[j])*this->U[j]) / (sqrt(p0[j] + y*ps[j]) + sqrt(q0[j] + y*qs[j]));
					}
				}

				//...Update Bl with BFGS...
				Vector<T> sl = alphay*dy;
				Vector<T> ql = (this->dWy(rs, ps, qs, _xk) - z) - dL;
				T psi = 1.0;
				if(sl*ql <= 0.2*sl*(Bl*sl)){
					psi = 0.8*sl*(Bl*sl) / (sl*(Bl*sl - ql));
				}
				Vector<T> qhat = psi*ql + (1.0 - psi)*(Bl*sl);
				Bl += -((Bl*sl)*((Bl*sl).Transpose())) / (sl*(Bl*sl)) + (qhat*(qhat.Transpose())) / (sl*qhat);
			}

			//...Update mu...
			std::cout << std::endl << y(0) << "\t" << z(0) << "\t" << mu;
			mu *= this->tau;
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
	Vector<T> MMA<T>::dWy(Vector<T> _rs, std::vector<Vector<T> > _ps, std::vector<Vector<T> > _qs, std::vector<T> _x) {
		Vector<T> vec = _rs;
		for(int j = 0; j < this->n; j++){
			vec += _ps[j] / (this->U[j] - _x[j]) + _qs[j] / (_x[j] - this->L[j]);
		}
		return -vec;
	}


	template<class T>
	T MMA<T>::LogBallier(T _mu, Vector<T> _y){
		T value = T();
		for(int i = 0; i < this->m; i++){
			value += -log(_y(i));
		}
		return (_mu)*value;
	}
}