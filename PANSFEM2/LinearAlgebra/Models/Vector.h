//*****************************************************************************
//Title		:PANSFEM2/LinearAlgebra/Models/Vector.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <vector>


namespace PANSFEM2 {
	//********************Vector Class********************
	template<class T>
	class Vector
	{
	public:
		Vector();
		~Vector();
		Vector(std::vector<T> _xs);
		Vector(int _size);


		const int DOX;
		std::vector<T> x;


		template<class F>
		friend std::ostream& operator << (std::ostream &_out, const Vector<F> &_p);
	};


	template<class T>
	inline Vector<T>::Vector() : DOX(0) {}


	template<class T>
	inline Vector<T>::~Vector() {}


	template<class T>
	inline Vector<T>::Vector(std::vector<T> _xs) : DOX(_xs.size()) {
		x = _xs;
	}


	template<class T>
	inline Vector<T>::Vector(int _size) : DOX(_size) {
		x = std::vector<T>(_size, T());
	}

	
	template<class F>
	std::ostream & operator<<(std::ostream & _out, const Vector<F>& _p) {
		for (auto xi : _p.x) {
			_out << xi << "\t";
		}
		//_out << std::endl;
		return _out;
	}


	//********************Vector Operations********************
	template<class T>
	T Triangle2DSpace(const Vector<T>& _p0, const Vector<T>& _p1, const Vector<T>& _p2) {
		return 0.5*((_p0.x[0] - _p2.x[0])*(_p1.x[1] - _p2.x[1]) - (_p2.x[1] - _p0.x[1])*(_p2.x[0] - _p1.x[0]));
	}
}