//*****************************************************************************
//Title		:src/FEM/Controller/ShapeFunction.h
//Author	:Tanabe Yuta
//Date		:2019/10/12
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************2NodesLine********************
	template<class T>
	class ShapeFunction2Line{
public:
		static const int d = 1;								//	Dimension on natural cordinate
		static const int n = 2;								//	Number of shape function
		static Vector<T> N(const std::vector<T>& _r);		//	Vector of shape function value
		static Matrix<T> dNdr(const std::vector<T>& _r);	//	Matrix of shape function derivative
	};


	template<class T>
	Vector<T> ShapeFunction2Line<T>::N(const std::vector<T>& _r){
		Vector<T> N = Vector<T>(n);
		N(0) = 0.5*(1 - _r[0]);
		N(1) = 0.5*(1 + _r[0]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction2Line<T>::dNdr(const std::vector<T>& _r){
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = -0.5;	dNdr(0, 1) = 0.5;
		return dNdr;
	}


	//********************3NodesLine********************
	template<class T>
	class ShapeFunction3Line{
public:
		static const int d = 1;
		static const int n = 3;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction3Line<T>::N(const std::vector<T>& _r){
		Vector<T> N = Vector<T>(n);
		N(0) = -0.5*(1.0 - _r[0])*_r[0];
		N(1) = 0.5*_r[0]*(1.0 + _r[0]);
		N(2) = (1.0 - _r[0])*(1.0 + _r[0]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction3Line<T>::dNdr(const std::vector<T>& _r){
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = -0.5*(1.0 - 2.0*_r[0]);	dNdr(0, 1) = 0.5*(1.0 + 2.0*_r[0]);		dNdr(0, 2) = -2.0*_r[0];
		return dNdr;
	}


	//********************3NodesTriangle********************
	template<class T>
	class ShapeFunction3Triangle {
public:
		static const int d = 2;
		static const int n = 3;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction3Triangle<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = _r[0];
		N(1) = _r[1];
		N(2) = 1.0 - _r[0] - _r[1];
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction3Triangle<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = 1.0;	dNdr(0, 1) = 0.0;	dNdr(0, 2) = -1.0;	
		dNdr(1, 0) = 0.0;	dNdr(1, 1) = 1.0;	dNdr(1, 2) = -1.0;	
		return dNdr;
	}


	//********************6NodesTriangle********************
	template<class T>
	class ShapeFunction6Triangle {
public:
		static const int d = 2;
		static const int n = 6;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction6Triangle<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = _r[0]*(2.0*_r[0] - 1.0);
		N(1) = _r[1]*(2.0*_r[1] - 1.0);
		N(2) = (1.0 - _r[0] - _r[1])*(1.0 - 2.0*_r[0] - 2.0*_r[1]);
		N(3) = 4.0*_r[0]*_r[1];
		N(4) = 4.0*_r[1]*(1.0 - _r[0] - _r[1]);
		N(5) = 4.0*(1.0 - _r[0] - _r[1])*_r[0];
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction6Triangle<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = 4.0*_r[0] - 1.0;	dNdr(0, 1) = 0.0;				dNdr(0, 2) = -3.0 + 4.0*_r[0] + 4.0*_r[1];	dNdr(0, 3) = 4.0*_r[1];	dNdr(0, 4) = -4.0*_r[1];					dNdr(0, 5) = 4.0*(1.0 - 2.0*_r[0] - _r[1]);
		dNdr(1, 0) = 0.0;				dNdr(1, 1) = 4.0*_r[1] - 1.0;	dNdr(1, 2) = -3.0 + 4.0*_r[0] + 4.0*_r[1];	dNdr(1, 3) = 4.0*_r[0];	dNdr(1, 4) = 4.0*(1.0 - _r[0] - 2.0*_r[1]);	dNdr(1, 5) = -4.0*_r[0];
		return dNdr;
	}


	//********************4NodesSquare********************
	template<class T>
	class ShapeFunction4Square {
public:
		static const int d = 2;
		static const int n = 4;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction4Square<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = 0.25*(1.0 - _r[0])*(1.0 - _r[1]);
		N(1) = 0.25*(1.0 + _r[0])*(1.0 - _r[1]);
		N(2) = 0.25*(1.0 + _r[0])*(1.0 + _r[1]);
		N(3) = 0.25*(1.0 - _r[0])*(1.0 + _r[1]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction4Square<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = -0.25*(1.0 - _r[1]);	dNdr(0, 1) = 0.25*(1.0 - _r[1]);	dNdr(0, 2) = 0.25*(1.0 - _r[1]);	dNdr(0, 3) = -0.25*(1.0 - _r[1]);	
		dNdr(1, 0) = -0.25*(1.0 - _r[0]);	dNdr(1, 1) = -0.25*(1.0 - _r[0]);	dNdr(1, 2) = 0.25*(1.0 - _r[0]);	dNdr(1, 3) = 0.25*(1.0 - _r[0]);
		return dNdr;
	}


	//********************8NodesSquare********************
	template<class T>
	class ShapeFunction8Square {
public:
		static const int d = 2;
		static const int n = 8;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction8Square<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = 0.25*(1.0 - _r[0])*(1.0 - _r[1])*(-_r[0] - _r[1] - 1.0);
		N(1) = 0.25*(1.0 + _r[0])*(1.0 - _r[1])*(_r[0] - _r[1] - 1.0);
		N(2) = 0.25*(1.0 + _r[0])*(1.0 + _r[1])*(_r[0] + _r[1] - 1.0);
		N(3) = 0.25*(1.0 - _r[0])*(1.0 + _r[1])*(-_r[0] + _r[1] - 1.0);
		N(4) = 0.5*(1.0 - _r[0])*(1.0 + _r[0])*(1.0 - _r[1]);
		N(5) = 0.5*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[1]);
		N(6) = 0.5*(1.0 + _r[0])*(1.0 - _r[0])*(1.0 + _r[1]);
		N(7) = 0.5*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[1]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction8Square<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = 0.25*(-(1.0 - _r[1])*(-_r[0] - _r[1] - 1.0) - (1.0 - _r[0])*(1.0 - _r[1]));
		dNdr(0, 1) = 0.25*((1.0 - _r[1])*(_r[0] - _r[1] - 1.0) + (1.0 + _r[0])*(1.0 - _r[1]));
		dNdr(0, 2) = 0.25*((1.0 + _r[1])*(_r[0] + _r[1] - 1.0) + (1.0 + _r[0])*(1.0 + _r[1]));
		dNdr(0, 3) = 0.25*(-(1.0 + _r[1])*(-_r[0] + _r[1] - 1.0) - (1.0 - _r[0])*(1.0 + _r[1]));
		dNdr(0, 4) = -_r[0] * (1.0 - _r[1]);
		dNdr(0, 5) = 0.5*(1.0 + _r[1])*(1.0 - _r[1]);
		dNdr(0, 6) = -_r[0] * (1.0 + _r[1]);
		dNdr(0, 7) = -0.5*(1.0 + _r[1])*(1.0 - _r[1]);

		dNdr(1, 0) = 0.25*(-(1.0 - _r[0])*(-_r[0] - _r[1] - 1.0) - (1.0 - _r[0])*(1.0 - _r[1]));
		dNdr(1, 1) = 0.25*(-(1.0 + _r[0])*(_r[0] - _r[1] - 1.0) - (1.0 + _r[0])*(1.0 - _r[1]));
		dNdr(1, 2) = 0.25*((1.0 + _r[0])*(_r[0] + _r[1] - 1.0) + (1.0 + _r[0])*(1.0 + _r[1]));
		dNdr(1, 3) = 0.25*((1.0 - _r[0])*(-_r[0] + _r[1] - 1.0) + (1.0 - _r[0])*(1.0 + _r[1]));
		dNdr(1, 4) = -0.5*(1.0 + _r[0])*(1.0 - _r[0]);
		dNdr(1, 5) = -_r[1] * (1.0 + _r[0]);
		dNdr(1, 6) = 0.5*(1.0 + _r[0])*(1.0 - _r[0]);
		dNdr(1, 7) = -_r[1] * (1.0 - _r[0]);		
		return dNdr;
	}


	//********************4NodesTetrahedron********************
	template<class T>
	class ShapeFunction4Tetrahedron {
public:
		static const int d = 3;
		static const int n = 4;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction4Tetrahedron<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = _r[0];
		N(1) = _r[1];
		N(2) = _r[2];
		N(3) = 1.0 - _r[0] - _r[1] - _r[2];
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction4Tetrahedron<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = 1.0;	dNdr(0, 1) = 0.0;	dNdr(0, 2) = 0.0;	dNdr(0, 3) = -1.0;
		dNdr(1, 0) = 0.0;	dNdr(1, 1) = 1.0;	dNdr(1, 2) = 0.0;	dNdr(1, 3) = -1.0;
		dNdr(2, 0) = 0.0;	dNdr(2, 1) = 0.0;	dNdr(2, 2) = 1.0;	dNdr(2, 3) = -1.0;
		return dNdr;
	}


	//********************8NodesCubic********************
	template<class T>
	class ShapeFunction8Cubic{
public:
		static const int d = 3;
		static const int n = 8;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction8Cubic<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = 0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 - _r[2]);
		N(1) = 0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 - _r[2]);
		N(2) = 0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[2]);
		N(3) = 0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[2]);
		N(4) = 0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[2]);
		N(5) = 0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 + _r[2]);
		N(6) = 0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 + _r[2]);
		N(7) = 0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[2]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction8Cubic<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = -0.125*(1.0 - _r[1])*(1.0 - _r[2]);	dNdr(0, 1) = 0.125*(1.0 - _r[1])*(1.0 - _r[2]);		dNdr(0, 2) = 0.125*(1.0 + _r[1])*(1.0 - _r[2]);		dNdr(0, 3) = -0.125*(1.0 + _r[1])*(1.0 - _r[2]);
		dNdr(0, 4) = -0.125*(1.0 - _r[1])*(1.0 + _r[2]);	dNdr(0, 5) = 0.125*(1.0 - _r[1])*(1.0 + _r[2]);		dNdr(0, 6) = 0.125*(1.0 + _r[1])*(1.0 + _r[2]);		dNdr(0, 7) = -0.125*(1.0 + _r[1])*(1.0 + _r[2]);

		dNdr(1, 0) = -0.125*(1.0 - _r[2])*(1.0 - _r[0]);	dNdr(1, 1) = -0.125*(1.0 - _r[2])*(1.0 + _r[0]);	dNdr(1, 2) = 0.125*(1.0 - _r[2])*(1.0 + _r[0]);		dNdr(1, 3) = 0.125*(1.0 - _r[2])*(1.0 - _r[0]);
		dNdr(1, 4) = -0.125*(1.0 + _r[2])*(1.0 - _r[0]);	dNdr(1, 5) = -0.125*(1.0 + _r[2])*(1.0 + _r[0]);	dNdr(1, 6) = 0.125*(1.0 + _r[2])*(1.0 + _r[0]);		dNdr(1, 7) = 0.125*(1.0 + _r[2])*(1.0 - _r[0]);

		dNdr(2, 0) = -0.125*(1.0 - _r[0])*(1.0 - _r[1]);	dNdr(2, 1) = -0.125*(1.0 + _r[0])*(1.0 - _r[1]);	dNdr(2, 2) = -0.125*(1.0 + _r[0])*(1.0 + _r[1]);	dNdr(2, 3) = -0.125*(1.0 - _r[0])*(1.0 + _r[1]);
		dNdr(2, 4) = 0.125*(1.0 - _r[0])*(1.0 - _r[1]);		dNdr(2, 5) = 0.125*(1.0 + _r[0])*(1.0 - _r[1]);		dNdr(2, 6) = 0.125*(1.0 + _r[0])*(1.0 + _r[1]);		dNdr(2, 7) = 0.125*(1.0 - _r[0])*(1.0 + _r[1]);
		return dNdr;
	}


	//********************20NodesIsoParametricElement********************
	template<class T>
	class ShapeFunction20Cubic{
public:
		static const int d = 3;
		static const int n = 20;
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction20Cubic<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(n);
		N(0) = -0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 - _r[2])*(2.0 + _r[0] + _r[1] + _r[2]);
		N(1) = -0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 - _r[2])*(2.0 - _r[0] + _r[1] + _r[2]);
		N(2) = -0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[2])*(2.0 - _r[0] - _r[1] + _r[2]);
		N(3) = -0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[2])*(2.0 + _r[0] - _r[1] + _r[2]);
		N(4) = -0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[2])*(2.0 + _r[0] + _r[1] - _r[2]);
		N(5) = -0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 + _r[2])*(2.0 - _r[0] + _r[1] - _r[2]);
		N(6) = -0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 + _r[2])*(2.0 - _r[0] - _r[1] - _r[2]);
		N(7) = -0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[2])*(2.0 + _r[0] - _r[1] - _r[2]);
		N(8) = 0.25*(1.0 + _r[0])*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 - _r[2]);
		N(9) = 0.25*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[1])*(1.0 - _r[2]);
		N(10) = 0.25*(1.0 + _r[0])*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[2]);
		N(11) = 0.25*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[1])*(1.0 - _r[2]);
		N(12) = 0.25*(1.0 + _r[0])*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[2]);
		N(13) = 0.25*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[1])*(1.0 + _r[2]);
		N(14) = 0.25*(1.0 + _r[0])*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[2]);
		N(15) = 0.25*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 - _r[1])*(1.0 + _r[2]);
		N(16) = 0.25*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[2])*(1.0 - _r[2]);
		N(17) = 0.25*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 + _r[2])*(1.0 - _r[2]);
		N(18) = 0.25*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 + _r[2])*(1.0 - _r[2]);
		N(19) = 0.25*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[2])*(1.0 - _r[2]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction20Cubic<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(d, n);
		dNdr(0, 0) = 0.125*(1.0 - _r[1])*(1.0 - _r[2])*(1.0 + 2.0*_r[0] + _r[1] + _r[2]);
		dNdr(0, 1) = -0.125*(1.0 - _r[1])*(1.0 - _r[2])*(1.0 - 2.0*_r[0] + _r[1] + _r[2]);
		dNdr(0, 2) = -0.125*(1.0 + _r[1])*(1.0 - _r[2])*(1.0 - 2.0*_r[0] - _r[1] + _r[2]);
		dNdr(0, 3) = 0.125*(1.0 + _r[1])*(1.0 - _r[2])*(1.0 + 2.0*_r[0] - _r[1] + _r[2]);
		dNdr(0, 4) = 0.125*(1.0 - _r[1])*(1.0 + _r[2])*(1.0 + 2.0*_r[0] + _r[1] - _r[2]);
		dNdr(0, 5) = -0.125*(1.0 - _r[1])*(1.0 + _r[2])*(1.0 - 2.0*_r[0] + _r[1] - _r[2]);
		dNdr(0, 6) = -0.125*(1.0 + _r[1])*(1.0 + _r[2])*(1.0 - 2.0*_r[0] - _r[1] - _r[2]);
		dNdr(0, 7) = 0.125*(1.0 + _r[1])*(1.0 + _r[2])*(1.0 + 2.0*_r[0] - _r[1] - _r[2]);
		dNdr(0, 8) = -0.5*_r[0]*(1.0 - _r[1])*(1.0 - _r[2]);
		dNdr(0, 9) = 0.25*(1.0 - _r[1] * _r[1])*(1.0 - _r[2]);
		dNdr(0, 10) = -0.5*_r[0]*(1.0 + _r[1])*(1.0 - _r[2]);
		dNdr(0, 11) = -0.25*(1.0 - _r[1] * _r[1])*(1.0 - _r[2]);
		dNdr(0, 12) = -0.5*_r[0]*(1.0 - _r[1])*(1.0 + _r[2]);
		dNdr(0, 13) = 0.25*(1.0 - _r[1] * _r[1])*(1.0 + _r[2]);
		dNdr(0, 14) = -0.5*_r[0]*(1.0 + _r[1])*(1.0 + _r[2]);
		dNdr(0, 15) = -0.25*(1.0 - _r[1] * _r[1])*(1.0 + _r[2]);
		dNdr(0, 16) = -0.25*(1.0 - _r[1])*(1.0 - _r[2] * _r[2]);
		dNdr(0, 17) = 0.25*(1.0 - _r[1])*(1.0 - _r[2] * _r[2]);
		dNdr(0, 18) = 0.25*(1.0 + _r[1])*(1.0 - _r[2] * _r[2]);
		dNdr(0, 19) = -0.25*(1.0 + _r[1])*(1.0 - _r[2] * _r[2]);

		dNdr(1, 0) = 0.125*(1.0 - _r[0])*(1.0 - _r[2])*(1.0 + _r[0] + 2.0*_r[1] + _r[2]);
		dNdr(1, 1) = 0.125*(1.0 + _r[0])*(1.0 - _r[2])*(1.0 - _r[0] + 2.0*_r[1] + _r[2]);
		dNdr(1, 2) = -0.125*(1.0 + _r[0])*(1.0 - _r[2])*(1.0 - _r[0] - 2.0*_r[1] + _r[2]);
		dNdr(1, 3) = -0.125*(1.0 - _r[0])*(1.0 - _r[2])*(1.0 + _r[0] - 2.0*_r[1] + _r[2]);
		dNdr(1, 4) = 0.125*(1.0 - _r[0])*(1.0 + _r[2])*(1.0 + _r[0] + 2.0*_r[1] - _r[2]);
		dNdr(1, 5) = 0.125*(1.0 + _r[0])*(1.0 + _r[2])*(1.0 - _r[0] + 2.0*_r[1] - _r[2]);
		dNdr(1, 6) = -0.125*(1.0 + _r[0])*(1.0 + _r[2])*(1.0 - _r[0] - 2.0*_r[1] - _r[2]);
		dNdr(1, 7) = -0.125*(1.0 - _r[0])*(1.0 + _r[2])*(1.0 + _r[0] - 2.0*_r[1] - _r[2]);
		dNdr(1, 8) = -0.25*(1.0 - _r[0] * _r[0])*(1.0 - _r[2]);
		dNdr(1, 9) = -0.5*(1.0 + _r[0])*_r[1]*(1.0 - _r[2]);
		dNdr(1, 10) = 0.25*(1.0 - _r[0] * _r[0])*(1.0 - _r[2]);
		dNdr(1, 11) = -0.5*(1.0 - _r[0])*_r[1]*(1.0 - _r[2]);
		dNdr(1, 12) = -0.25*(1.0 - _r[0] * _r[0])*(1.0 + _r[2]);
		dNdr(1, 13) = -0.5*(1.0 + _r[0])*_r[1]*(1.0 + _r[2]);
		dNdr(1, 14) = 0.25*(1.0 - _r[0] * _r[0])*(1.0 + _r[2]);
		dNdr(1, 15) = -0.5*(1.0 - _r[0])*_r[1]*(1.0 + _r[2]);
		dNdr(1, 16) = -0.25*(1.0 - _r[0])*(1.0 - _r[2] * _r[2]);
		dNdr(1, 17) = -0.25*(1.0 + _r[0])*(1.0 - _r[2] * _r[2]);
		dNdr(1, 18) = 0.25*(1.0 + _r[0])*(1.0 - _r[2] * _r[2]);
		dNdr(1, 19) = 0.25*(1.0 - _r[0])*(1.0 - _r[2] * _r[2]);

		dNdr(2, 0) = 0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[0] + _r[1] + 2.0*_r[2]);
		dNdr(2, 1) = 0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 - _r[0] + _r[1] + 2.0*_r[2]);
		dNdr(2, 2) = 0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[0] - _r[1] + 2.0*_r[2]);
		dNdr(2, 3) = 0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[0] - _r[1] + 2.0*_r[2]);
		dNdr(2, 4) = -0.125*(1.0 - _r[0])*(1.0 - _r[1])*(1.0 + _r[0] + _r[1] - 2.0*_r[2]);
		dNdr(2, 5) = -0.125*(1.0 + _r[0])*(1.0 - _r[1])*(1.0 - _r[0] + _r[1] - 2.0*_r[2]);
		dNdr(2, 6) = -0.125*(1.0 + _r[0])*(1.0 + _r[1])*(1.0 - _r[0] - _r[1] - 2.0*_r[2]);
		dNdr(2, 7) = -0.125*(1.0 - _r[0])*(1.0 + _r[1])*(1.0 + _r[0] - _r[1] - 2.0*_r[2]);
		dNdr(2, 8) = -0.25*(1.0 - _r[0] * _r[0])*(1.0 - _r[1]);
		dNdr(2, 9) = -0.25*(1.0 + _r[0])*(1.0 - _r[1] * _r[1]);
		dNdr(2, 10) = -0.25*(1.0 - _r[0] * _r[0])*(1.0 + _r[1]);
		dNdr(2, 11) = -0.25*(1.0 - _r[0])*(1.0 - _r[1] * _r[1]);
		dNdr(2, 12) = 0.25*(1.0 - _r[0] * _r[0])*(1.0 - _r[1]);
		dNdr(2, 13) = 0.25*(1.0 + _r[0])*(1.0 - _r[1] * _r[1]);
		dNdr(2, 14) = 0.25*(1.0 - _r[0] * _r[0])*(1.0 + _r[1]);
		dNdr(2, 15) = 0.25*(1.0 - _r[0])*(1.0 - _r[1] * _r[1]);
		dNdr(2, 16) = -0.5*(1.0 - _r[0])*(1.0 - _r[1])*_r[2];
		dNdr(2, 17) = -0.5*(1.0 + _r[0])*(1.0 - _r[1])*_r[2];
		dNdr(2, 18) = -0.5*(1.0 + _r[0])*(1.0 + _r[1])*_r[2];
		dNdr(2, 19) = -0.5*(1.0 - _r[0])*(1.0 + _r[1])*_r[2];		
		return dNdr;
	}
}