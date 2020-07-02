//*********************************************************
//  Title       :   src/FEM/Equation/ReactionDiffusion.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/05/04
//  Copyright   :   (C)2020 TanabeYuta
//*********************************************************


#pragma once
#include <vector>
#include <cassert>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Reaction-Diffusion consistent mass matrix********************
    template<class T, template<class>class SF, template<class>class IC>
    void ReactionDiffusionConsistentMass(Matrix<T>& _Me, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x) {
        assert(_doulist.size() == 1);

		_Me = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

        Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);   X(i, 1) = _x[_element[i]](1);
		}

        for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
            _Me += N*N.Transpose()*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }
    }


    //********************Reaction-Diffusion diffusion matrix********************
    template<class T, template<class>class SF, template<class>class IC>
    void ReactionDiffusionStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _D) {
        assert(_doulist.size() == 1);

        _Ke = Matrix<T>(_element.size(), _element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

        Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);   X(i, 1) = _x[_element[i]](1);
		}

        for (int g = 0; g < IC<T>::N; g++) {
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
            _Ke += _D*dNdX.Transpose()*dNdX*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }
    }


    //********************Reaction-Diffusion reaction vector********************
    template<class T, template<class>class SF, template<class>class IC, class F>
    void ReactionDiffusionReaction(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, std::vector<Vector<T> >& _u, F _f) {
        assert(_doulist.size() == 1);

        _Fe = Vector<T>(_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

        Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);   X(i, 1) = _x[_element[i]](1);
		}

        Vector<T> U = Vector<T>(_element.size());
        for(int i = 0; i < _element.size(); i++) {
            U(i) = _u[_element[i]](0);
        }

        for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

            T u = N*U;
            Vector<T> dudX = dNdX*U;
           
            _Fe += N*_f(u, dudX)*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }
    }





    //********************Make LevelSet method equation********************
    template<class T, template<class>class SF, template<class>class IC>
    void LevelSet(Matrix<T>& _Te, Vector<T>& _Ye, std::vector<std::vector<std::pair<int, int> > >& _nodetoelement, const std::vector<int>& _element, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _dt, T _tau, std::vector<T>& _F, std::vector<Vector<T> >& _phi) {
        assert(_doulist.size() == 1);

		_Te = Matrix<T>(_element.size(), _element.size());
        _Ye = Vector<T>(_element.size());
		_nodetoelement = std::vector<std::vector<std::pair<int, int> > >(_element.size(), std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < _element.size(); i++) {
			_nodetoelement[i][0] = std::make_pair(_doulist[0], i);
		}

        Matrix<T> X = Matrix<T>(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);
			X(i, 1) = _x[_element[i]](1);
		}

        Vector<T> Phi = Vector<T>(_element.size());
        for(int i = 0; i < _element.size(); i++) {
            Phi(i) = _phi[_element[i]](0);
        }

        Vector<T> F = Vector<T>(_element.size());
        for(int i = 0; i < _element.size(); i++) {
            F(i) = _F[_element[i]];
        } 

        for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
            T phi = N*Phi;
            T f = N*F; 
            _Te += (N*N.Transpose()/_dt + _tau*dNdX.Transpose()*dNdX)*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
            _Ye += (f + phi/_dt)*N*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }
    }
}