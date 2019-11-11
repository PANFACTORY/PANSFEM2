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
	//********************3NodesTriangle********************
	template<class T>
	class ShapeFunction3Triangle {
public:
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction3Triangle<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(3);
		N(0) = _r[0];
		N(1) = _r[1];
		N(2) = 1.0 - _r[0] - _r[1];
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction3Triangle<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(2, 3);
		dNdr(0, 0) = 1.0;	dNdr(0, 1) = 0.0;	dNdr(0, 2) = -1.0;	
		dNdr(1, 0) = 0.0;	dNdr(1, 1) = 1.0;	dNdr(1, 2) = -1.0;	
		return dNdr;
	}


	//********************4NodesSquare********************
	template<class T>
	class ShapeFunction4Square {
public:
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction4Square<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(4);
		N(0) = 0.25*(1.0 - _r[0])*(1.0 - _r[1]);
		N(1) = 0.25*(1.0 + _r[0])*(1.0 - _r[1]);
		N(2) = 0.25*(1.0 + _r[0])*(1.0 + _r[1]);
		N(3) = 0.25*(1.0 - _r[0])*(1.0 + _r[1]);
		return N;
	}


	template<class T>
	Matrix<T> ShapeFunction4Square<T>::dNdr(const std::vector<T>& _r) {
		Matrix<T> dNdr = Matrix<T>(2, 4);
		dNdr(0, 0) = -0.25*(1.0 - _r[1]);	dNdr(0, 1) = 0.25*(1.0 - _r[1]);	dNdr(0, 2) = 0.25*(1.0 - _r[1]);	dNdr(0, 3) = -0.25*(1.0 - _r[1]);	
		dNdr(1, 0) = -0.25*(1.0 - _r[0]);	dNdr(1, 1) = -0.25*(1.0 - _r[0]);	dNdr(1, 2) = 0.25*(1.0 - _r[0]);	dNdr(1, 3) = 0.25*(1.0 - _r[0]);
		return dNdr;
	}


	//********************8NodesSquare********************
	template<class T>
	class ShapeFunction8Square {
public:
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction8Square<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(8);
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
		Matrix<T> dNdr = Matrix<T>(2, 8);
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


	//********************8NodesCubic********************
	template<class T>
	class ShapeFunction8Cubic{
public:
		static Vector<T> N(const std::vector<T>& _r);
		static Matrix<T> dNdr(const std::vector<T>& _r);
	};


	template<class T>
	Vector<T> ShapeFunction8Cubic<T>::N(const std::vector<T>& _r) {
		Vector<T> N = Vector<T>(8);
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
		Matrix<T> dNdr = Matrix<T>(3, 8);
		dNdr(0, 0) = -0.125*(1.0 - _r[1])*(1.0 - _r[2]);	dNdr(0, 1) = 0.125*(1.0 - _r[1])*(1.0 - _r[2]);		dNdr(0, 2) = 0.125*(1.0 + _r[1])*(1.0 - _r[2]);		dNdr(0, 3) = -0.125*(1.0 + _r[1])*(1.0 - _r[2]);
		dNdr(0, 4) = -0.125*(1.0 - _r[1])*(1.0 + _r[2]);	dNdr(0, 5) = 0.125*(1.0 - _r[1])*(1.0 + _r[2]);		dNdr(0, 6) = 0.125*(1.0 + _r[1])*(1.0 + _r[2]);		dNdr(0, 7) = -0.125*(1.0 + _r[1])*(1.0 + _r[2]);

		dNdr(1, 0) = -0.125*(1.0 - _r[2])*(1.0 - _r[0]);	dNdr(1, 1) = -0.125*(1.0 - _r[2])*(1.0 + _r[0]);	dNdr(1, 2) = 0.125*(1.0 - _r[2])*(1.0 + _r[0]);		dNdr(1, 3) = 0.125*(1.0 - _r[2])*(1.0 - _r[0]);
		dNdr(1, 4) = -0.125*(1.0 + _r[2])*(1.0 - _r[0]);	dNdr(1, 5) = -0.125*(1.0 + _r[2])*(1.0 + _r[0]);	dNdr(1, 6) = 0.125*(1.0 + _r[2])*(1.0 + _r[0]);		dNdr(1, 7) = 0.125*(1.0 + _r[2])*(1.0 - _r[0]);

		dNdr(2, 0) = -0.125*(1.0 - _r[0])*(1.0 - _r[1]);	dNdr(2, 1) = -0.125*(1.0 + _r[0])*(1.0 - _r[1]);	dNdr(2, 2) = -0.125*(1.0 + _r[0])*(1.0 + _r[1]);	dNdr(2, 3) = -0.125*(1.0 - _r[0])*(1.0 + _r[1]);
		dNdr(2, 4) = 0.125*(1.0 - _r[0])*(1.0 - _r[1]);		dNdr(2, 5) = 0.125*(1.0 + _r[0])*(1.0 - _r[1]);		dNdr(2, 6) = 0.125*(1.0 + _r[0])*(1.0 + _r[1]);		dNdr(2, 7) = 0.125*(1.0 - _r[0])*(1.0 + _r[1]);
		return dNdr;
	}


	//********************20NodesIsoParametricElement********************
	
}