//*****************************************************************************
//  Title       :src/LinearAlgebra/Solvers/Lanczos.h
//  Author      :Tanabe Yuta
//  Date        :2019/11/22
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <numeric>


#include "../Models/CSR.h"
#include "CG.h"


//********************{x}={x}/a********************
template<class T>
void xexda(std::vector<T>& _x, T _a) {
    for(auto& xi : _x) {
        xi /= _a;
    }
}


//********************{x}={y}/a********************
template<class T>
void xeyda(std::vector<T>& _x, const std::vector<T>& _y, T _a) {
    auto yi = _y.begin();
    for(auto& xi : _x) {
        xi = (*yi)/_a;
        ++yi;
    }
}


//********************Bisection method********************
template<class T>
T BisectionMethod(const std::vector<T>& _alpha, const std::vector<T>& _beta, int _mode){
    int m = _alpha.size();
    
    T b = fabs(_alpha[0]) + fabs(_beta[0]);
    for(int i = 1; i < m; i++){
        if(b < fabs(_beta[i - 1]) + fabs(_alpha[i]) + fabs(_beta[i])){
            b = fabs(_beta[i - 1]) + fabs(_alpha[i]) + fabs(_beta[i]);
        }
    }

    T a1 = -b, a2 = b, lambda;
    while((a2 - a1) > 1.0e-13*b){
        //----------Update lambda----------
        lambda = 0.5*(a1 + a2);

        //----------Count sign change----------
        T pkm1 = 1.0, pk = _alpha[0] - lambda, sign = 1.0;
        int Nlambda = 0;
        if(sign*pk < T()){
            Nlambda++;
            sign *= -1.0;
        }
        for(int k = 0; k < m - 1; k++){
            T pkp1 = (_alpha[k + 1] - lambda)*pk - pow(_beta[k], 2.0)*pkm1;
            if(sign*pkp1 < T()){
                Nlambda++;
                sign *= -1.0;
            }
            pkm1 = pk;
            pk = pkp1;
        }

        //----------Update a1 and a2----------
        if(Nlambda <= _mode){
            a1 = lambda;
        } else {
            a2 = lambda;
        }
    }

    return lambda;
}


//********************Inverse Power method********************
template<class T>
std::vector<T> InversePowerMethod(const std::vector<T>& _alpha, const std::vector<T>& _beta, T _lambda){
    int m = _alpha.size();
    std::vector<T> yk = std::vector<T>(m, T());
    yk[0] = 1.0;

    std::vector<T> p = std::vector<T>(m);
    std::vector<T> q = std::vector<T>(m);
    for(int k = 0; k < 10000; k++){
        //.....Solve [T - lambda*E]{ykp1}={yk} with Thomas method.....
        p[0] = -_beta[0]/(_alpha[0] - _lambda);
        q[0] = yk[0]/(_alpha[0] - _lambda);
        for(int i = 1; i < m; i++){
            p[i] = -_beta[i]/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
            q[i] = (yk[i] - _beta[i - 1]*q[i - 1])/((_alpha[i] - _lambda) + _beta[i - 1]*p[i - 1]);
        }

        yk[m - 1] = q[m - 1];
        for(int i = m - 2; i >= 0; i--){
            yk[i] = p[i]*yk[i + 1] + q[i];
        }

        //.....Normalize yk.....
        T ykNorm = sqrt(std::inner_product(yk.begin(), yk.end(), yk.begin(), T()));
        std::transform(yk.begin(), yk.end(), yk.begin(), [=](T _yki) { return _yki/ykNorm; });
    }
    
    return yk;
}


//********************Reconvert vector********************
template<class T>
std::vector<T> ReconvertVector(const std::vector<T>& _y, std::vector<std::vector<T> >& _q){
    int n = _q[0].size();
    int m = _y.size();
    std::vector<T> x = std::vector<T>(n, T());
    for(int k = 0; k < m; k++){
        for(int i = 0; i < n; i++){
            x[i] += _q[k][i]*_y[k];
        } 
    }
    return x;
}


//********************Lanczos process********************
template<class T>
void Lanczos(CSR<T>& _A, std::vector<T>& _eigenvalues, std::vector<std::vector<T> >& _eigenvectors, int _m){
    //----------Initialize----------
    _eigenvalues = std::vector<T>(_m);                                                                  //  Eigenvalues
    _eigenvectors = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS));                          //  Eigenvectors          
    std::vector<T> alpha = std::vector<T>(_m);                                                          //  Values of diagonal
    std::vector<T> beta = std::vector<T>(_m);                                                           //  Values of side of diagonal
    std::vector<std::vector<T> > q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));    //  Orthogonal vectors
    q[0] = GetDiagonal(_A);
    xexda(q[0], sqrt(std::inner_product(q[0].begin(), q[0].end(), q[0].begin(), T())));

    //----------Lanczos process----------
    for(int k = 0; k < _m; k++){
        std::vector<T> p = _A*q[k];
        if( k != 0){
            xexpay(p, -beta[k - 1], q[k - 1]);
        }
        alpha[k] = std::inner_product(q[k].begin(), q[k].end(), p.begin(), T());
        xexpay(p, -alpha[k], q[k]);
        beta[k] = sqrt(std::inner_product(p.begin(), p.end(), p.begin(), T()));
        if(k != _m - 1){
            xeyda(q[k + 1], p, beta[k]);
        }
    }

    //----------Get eigenvalues and eigenvectors----------
    for(int k = 0; k < _m; k++){
		_eigenvalues[k] = BisectionMethod(alpha, beta, k);		
		_eigenvectors[k] = ReconvertVector(InversePowerMethod(alpha, beta, _eigenvalues[k]), q);
	}	
}


//********************Restart Lanczos process********************
template<class T>
void RestartLanczos(CSR<T>& _A, std::vector<T>& _eigenvalues, std::vector<std::vector<T> >& _eigenvectors, int _m){
    //----------Initialize----------
    _eigenvalues = std::vector<T>(_m);                                                                  //  Eigenvalues
    _eigenvectors = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS));                          //  Eigenvectors          
    std::vector<T> alpha = std::vector<T>(_m);                                                          //  Values of diagonal
    std::vector<T> beta = std::vector<T>(_m);                                                           //  Values of side of diagonal
    std::vector<std::vector<T> > q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));    //  Orthogonal vectors
    q[0] = GetDiagonal(_A);
    xexda(q[0], sqrt(std::inner_product(q[0].begin(), q[0].end(), q[0].begin(), T())));

    //----------Lanczos process----------
    for(int k = 0; k < _m; k++){
        std::vector<T> p = _A*q[k];
        for(int i = 0; i < k - 2; i++) {
            xexpay(p, -std::inner_product(p.begin(), p.end(), q[i].begin(), T()), q[i]);
        }
        if( k != 0){
            xexpay(p, -beta[k - 1], q[k - 1]);
        }
        alpha[k] = std::inner_product(q[k].begin(), q[k].end(), p.begin(), T());
        xexpay(p, -alpha[k], q[k]);
        beta[k] = sqrt(std::inner_product(p.begin(), p.end(), p.begin(), T()));
        if(k != _m - 1){
            xeyda(q[k + 1], p, beta[k]);
        }
    }

    //----------Get eigenvalues and eigenvectors----------
    for(int k = 0; k < _m; k++){
		_eigenvalues[k] = BisectionMethod(alpha, beta, k);		
		_eigenvectors[k] = ReconvertVector(InversePowerMethod(alpha, beta, _eigenvalues[k]), q);
	}	
}


//********************Shifted-Invert Lanczos process********************
template<class T>
void ShiftedInvertLanczos(CSR<T>& _A, std::vector<T>& _eigenvalues, std::vector<std::vector<T> >& _eigenvectors, int _m, T _sigma){  
    //----------Initialize----------
    LILCSR<T> I = LILCSR<T>(_A.ROWS, _A.COLS);
    for(int i = 0; i < _A.ROWS; i++) {
        I.set(i, i, _sigma);
    }
    CSR<T> A = _A - CSR<T>(I);
    _eigenvalues = std::vector<T>(_m);                                                                  //  Eigenvalues
    _eigenvectors = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS));                          //  Eigenvectors          
    std::vector<T> alpha = std::vector<T>(_m);                                                          //  Values of diagonal
    std::vector<T> beta = std::vector<T>(_m);                                                           //  Values of side of diagonal
    std::vector<std::vector<T> > q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));    //  Orthogonal vectors
    q[0] = GetDiagonal(A);
    xexda(q[0], sqrt(std::inner_product(q[0].begin(), q[0].end(), q[0].begin(), T())));
    int itrmax = std::max(_A.ROWS, 1000);

    //----------Lanczos process----------
    for(int k = 0; k < _m; k++){
        std::vector<T> p = ScalingCG(A, q[k], itrmax, 1.0e-10);
        if( k != 0){
            xexpay(p, -beta[k - 1], q[k - 1]);
        }
        alpha[k] = std::inner_product(q[k].begin(), q[k].end(), p.begin(), T());
        xexpay(p, -alpha[k], q[k]);
        beta[k] = sqrt(std::inner_product(p.begin(), p.end(), p.begin(), T()));
        if(k != _m - 1){
            xeyda(q[k + 1], p, beta[k]);
        }
    }

    //----------Get eigenvalues and eigenvectors----------
    for(int k = 0; k < _m; k++){
        T eigenvalue = BisectionMethod(alpha, beta, _m - k - 1);
		_eigenvalues[k] = _sigma + 1.0/eigenvalue;		
		_eigenvectors[k] = ReconvertVector(InversePowerMethod(alpha, beta, eigenvalue), q);
	}	
}


//********************Shifted-Invert Lanczos process for General eigenvalue problem********************
template<class T>
void GeneralShiftedInvertLanczos(CSR<T>& _A, CSR<T>& _B, std::vector<T>& _eigenvalues, std::vector<std::vector<T> >& _eigenvectors, int _m, T _sigma){
    //----------Initialize----------
    CSR<T> A = _A - _B*_sigma;
    _eigenvalues = std::vector<T>(_m);                                                                  //  Eigenvalues
    _eigenvectors = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS));                          //  Eigenvectors          
    std::vector<T> alpha = std::vector<T>(_m);                                                          //  Values of diagonal
    std::vector<T> beta = std::vector<T>(_m);                                                           //  Values of side of diagonal
    std::vector<std::vector<T> > q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));    //  Orthogonal vectors
    q[0] = GetDiagonal(A);
    xexda(q[0], sqrt(std::inner_product(q[0].begin(), q[0].end(), q[0].begin(), T())));
    std::vector<T> p = _B*q[0];
    int itrmax = std::max(_A.ROWS, 10000);

    //----------Lanczos process----------
    for(int k = 0; k < _m; k++){
        std::vector<T> s = ScalingCG(A, p, itrmax, 1.0e-10);
        if(k != 0){
            xexpay(s, -beta[k - 1], q[k - 1]);
        }
        alpha[k] = std::inner_product(p.begin(), p.end(), s.begin(), T());
        xexpay(s, -alpha[k], q[k]);
        std::vector<T> r = _B*s;
        beta[k] = sqrt(std::inner_product(r.begin(), r.end(), s.begin(), T()));
        xeyda(p, r, beta[k]);
        if(k != _m - 1){
            xeyda(q[k + 1], s, beta[k]);
        }       
    }
     
    //----------Get eigenvalues and eigenvectors----------
    for(int k = 0; k < _m; k++){
        T eigenvalue = BisectionMethod(alpha, beta, _m - k - 1);
		_eigenvalues[k] = _sigma + 1.0/eigenvalue;		
		_eigenvectors[k] = ReconvertVector(InversePowerMethod(alpha, beta, eigenvalue), q);
	}	
}


//********************Restart Invert Lanczos process for General eigenvalue problem********************
template<class T>
void GeneralRestartShiftedInvertLanczos(CSR<T>& _A, CSR<T>& _B, std::vector<T>& _eigenvalues, std::vector<std::vector<T> >& _eigenvectors, int _m, T _sigma){
    //----------Initialize----------
    CSR<T> A = _A - _B*_sigma;
    _eigenvalues = std::vector<T>(_m);                                                                  //  Eigenvalues
    _eigenvectors = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS));                          //  Eigenvectors          
    std::vector<T> alpha = std::vector<T>(_m);                                                          //  Values of diagonal
    std::vector<T> beta = std::vector<T>(_m);                                                           //  Values of side of diagonal
    std::vector<std::vector<T> > q = std::vector<std::vector<T> >(_m, std::vector<T>(_A.ROWS, T()));    //  Orthogonal vectors
    q[0] = GetDiagonal(A);
    xexda(q[0], sqrt(std::inner_product(q[0].begin(), q[0].end(), q[0].begin(), T())));
    std::vector<T> p = _B*q[0];
    int itrmax = std::max(_A.ROWS, 10000);

    //----------Lanczos process----------
    for(int k = 0; k < _m; k++){
        std::vector<T> s = ScalingCG(A, p, itrmax, 1.0e-10);
        for(int i = 0; i < k - 2; i++) {
            std::vector<T> Bq = _B*q[i];
            xexpay(s, -std::inner_product(s.begin(), s.end(), Bq.begin(), T()), Bq);
        }
        if(k != 0){
            xexpay(s, -beta[k - 1], q[k - 1]);
        }
        alpha[k] = std::inner_product(p.begin(), p.end(), s.begin(), T());
        xexpay(s, -alpha[k], q[k]);
        std::vector<T> r = _B*s;
        beta[k] = sqrt(std::inner_product(r.begin(), r.end(), s.begin(), T()));
        xeyda(p, r, beta[k]);
        if(k != _m - 1){
            xeyda(q[k + 1], s, beta[k]);
        }       
    }
     
    //----------Get eigenvalues and eigenvectors----------
    for(int k = 0; k < _m; k++){
        T eigenvalue = BisectionMethod(alpha, beta, _m - k - 1);
		_eigenvalues[k] = _sigma + 1.0/eigenvalue;		
		_eigenvectors[k] = ReconvertVector(InversePowerMethod(alpha, beta, eigenvalue), q);
	}	
}