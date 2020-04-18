//*****************************************************************************
//  Title		:   src/FEM/Equation/Beam.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/02/19
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Linear Isotropic Elastic Beam********************
	template<class T, template<class>class SF, template<class>class IC0, template<class>class IC12>
	void BeamLinearIsotropicElastic(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _v1, std::vector<Vector<T> >& _v2, T _E, T _V, T _a, T _b) {
        assert(_doulist.size() == 6);
        
        _Ke = Matrix<T>(6*_element.size(), 6*_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(6));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], 6*i);
			_nodetoelement[i][1] = std::make_pair(_doulist[1], 6*i + 1);
			_nodetoelement[i][2] = std::make_pair(_doulist[2], 6*i + 2);
            _nodetoelement[i][3] = std::make_pair(_doulist[3], 6*i + 3);
			_nodetoelement[i][4] = std::make_pair(_doulist[4], 6*i + 4);
            _nodetoelement[i][5] = std::make_pair(_doulist[5], 6*i + 5);
		}
		
       Matrix<T> X = Matrix<T>(_element.size(), 3);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);	X(i, 1) = _x[_element[i]](1);	X(i, 2) = _x[_element[i]](2);
		}

        Matrix<T> v1 = Matrix<T>(0, 3);
        Matrix<T> v2 = Matrix<T>(0, 3);
        for(auto i : _element){
            v1 = v1.Vstack(_v1[i].Normal().Transpose());
            v2 = v2.Vstack(_v2[i].Normal().Transpose());
        }

        T k = 1.2;
        Matrix<T> D0 = Matrix<T>(3, 3);
        D0(0, 0) = _E;  D0(0, 1) = T();                 D0(0, 2) = T();
        D0(1, 0) = T(); D0(1, 1) = 0.5*_E/(1.0 + _V)/k; D0(1, 2) = T();
        D0(2, 0) = T(); D0(2, 1) = T();                 D0(2, 2) = 0.5*_E/(1.0 + _V)/k;

        for(int g = 0; g < IC0<T>::N; g++){
            for(int h = 0; h < IC12<T>::N; h++){
                Vector<T> N = SF<T>::N(IC0<T>::Points[g]);
                Matrix<T> dNdr = SF<T>::dNdr(IC0<T>::Points[g]);
                
                Vector<T> itazeta = IC12<T>::Points[h];
                Matrix<T> J = (dNdr*(X + 0.5*_a*itazeta(0)*v1 + 0.5*_b*itazeta(1)*v2)).Vstack((0.5*_a*N.Transpose()*v1).Vstack(0.5*_b*N.Transpose()*v2));
                T detJ = J.Determinant();
                Matrix<T> invJ = J.Inverse();
                Matrix<T> dNdx = invJ.Block(0, 0, 3, 1)*dNdr;
                Matrix<T> dNydx = dNdx*itazeta(0) + invJ.Block(0, 1, 3, 1)*N.Transpose();
                Matrix<T> dNzdx = dNdx*itazeta(1) + invJ.Block(0, 2, 3, 1)*N.Transpose();
               
                Matrix<T> B = Matrix<T>(6, 6*_element.size());
                for(int i = 0; i < _element.size(); i++){
                    B(0, 6*i) = dNdx(0, i); B(0, 6*i + 1) = T();        B(0, 6*i + 2) = T();        B(0, 6*i + 3) = T();                                                                                                                    B(0, 6*i + 4) =  0.5*_a*dNydx(0, i)*v1(i, 2) + 0.5*_b*dNzdx(0, i)*v2(i, 2);                                                             B(0, 6*i + 5) = -0.5*_a*dNydx(0, i)*v1(i, 1) - 0.5*_b*dNzdx(0, i)*v2(i, 1);
                    B(1, 6*i) = T();        B(1, 6*i + 1) = dNdx(1, i); B(1, 6*i + 2) = T();        B(1, 6*i + 3) = -0.5*_a*dNydx(1, i)*v1(i, 2) - 0.5*_b*dNzdx(1, i)*v2(i, 2);                                                             B(1, 6*i + 4) = T();                                                                                                                    B(1, 6*i + 5) =  0.5*_a*dNydx(1, i)*v1(i, 0) + 0.5*_b*dNzdx(1, i)*v2(i, 0);
                    B(2, 6*i) = T();        B(2, 6*i + 1) = T();        B(2, 6*i + 2) = dNdx(2, i); B(2, 6*i + 3) =  0.5*_a*dNydx(2, i)*v1(i, 1) + 0.5*_b*dNzdx(2, i)*v2(i, 1);                                                             B(2, 6*i + 4) = -0.5*_a*dNydx(2, i)*v1(i, 0) - 0.5*_b*dNzdx(2, i)*v2(i, 0);                                                             B(2, 6*i + 5) = T();
                    B(3, 6*i) = dNdx(1, i); B(3, 6*i + 1) = dNdx(0, i); B(3, 6*i + 2) = T();        B(3, 6*i + 3) = -0.5*_a*dNydx(0, i)*v1(i, 2) - 0.5*_b*dNzdx(0, i)*v2(i, 2);                                                             B(3, 6*i + 4) =  0.5*_a*dNydx(1, i)*v1(i, 2) + 0.5*_b*dNzdx(1, i)*v2(i, 2);                                                             B(3, 6*i + 5) = -0.5*_a*dNydx(1, i)*v1(i, 1) - 0.5*_b*dNzdx(1, i)*v2(i, 1) + 0.5*_a*dNydx(0, i)*v1(i, 0) + 0.5*_b*dNzdx(0, i)*v2(i, 0);
                    B(4, 6*i) = dNdx(2, i); B(4, 6*i + 1) = T();        B(4, 6*i + 2) = dNdx(0, i); B(4, 6*i + 3) =  0.5*_a*dNydx(0, i)*v1(i, 1) + 0.5*_b*dNzdx(0, i)*v2(i, 1);                                                             B(4, 6*i + 4) = -0.5*_a*dNydx(0, i)*v1(i, 0) - 0.5*_b*dNzdx(0, i)*v2(i, 0) + 0.5*_a*dNydx(2, i)*v1(i, 2) + 0.5*_b*dNzdx(2, i)*v2(i, 2); B(4, 6*i + 5) = -0.5*_a*dNydx(2, i)*v1(i, 1) - 0.5*_b*dNzdx(2, i)*v2(i, 1);
                    B(5, 6*i) = T();        B(5, 6*i + 1) = dNdx(2, i); B(5, 6*i + 2) = dNdx(1, i); B(5, 6*i + 3) = -0.5*_a*dNydx(2, i)*v1(i, 2) - 0.5*_b*dNzdx(2, i)*v2(i, 2) + 0.5*_a*dNydx(1, i)*v1(i, 1) + 0.5*_b*dNzdx(1, i)*v2(i, 1); B(5, 6*i + 4) = -0.5*_a*dNydx(1, i)*v1(i, 0) - 0.5*_b*dNzdx(1, i)*v2(i, 0);                                                             B(5, 6*i + 5) =  0.5*_a*dNydx(2, i)*v1(i, 0) + 0.5*_b*dNzdx(2, i)*v2(i, 0);
                }

                Vector<T> P1 = (v1.Transpose()*N).Normal();
                Vector<T> P2 = (v2.Transpose()*N).Normal();
                Vector<T> P0 = VectorProduct(P1, P2).Normal();
                Matrix<T> P = Matrix<T>(3, 6);
                P(0, 0) = P0(0)*P0(0);      P(0, 1) = P0(1)*P0(1);      P(0, 2) = P0(2)*P0(2);      P(0, 3) = P0(0)*P0(1);                  P(0, 4) = P0(0)*P0(2);                  P(0, 5) = P0(1)*P0(2);
                P(1, 0) = 2.0*P0(0)*P1(0);  P(1, 1) = 2.0*P0(1)*P1(1);  P(1, 2) = 2.0*P0(2)*P1(2);  P(1, 3) = P0(0)*P1(1) + P0(1)*P1(0);    P(1, 4) = P0(0)*P1(2) + P0(2)*P1(0);    P(1, 5) = P0(1)*P1(2) + P0(2)*P1(1);
                P(2, 0) = 2.0*P0(0)*P2(0);  P(2, 1) = 2.0*P0(1)*P2(1);  P(2, 2) = 2.0*P0(2)*P2(2);  P(2, 3) = P0(0)*P2(1) + P0(1)*P2(0);    P(2, 4) = P0(0)*P2(2) + P0(2)*P2(0);    P(2, 5) = P0(1)*P2(2) + P0(2)*P2(1);
                Matrix<T> D = P.Transpose()*D0*P;

                _Ke += B.Transpose()*D*B*detJ*IC0<T>::Weights[g][0]*IC12<T>::Weights[h][0]*IC12<T>::Weights[h][1];
            }
        } 
    }
}