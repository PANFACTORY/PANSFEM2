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
	template<class T>
	class Vector
	{
	public:
		Vector();
		~Vector();
		template<class... Ts>
		Vector(Ts... _xs);
		Vector(std::vector<T> _xs);


		const int DOX;
		std::vector<T> x;


		template<class F>
		friend std::ostream& operator << (std::ostream &_out, const Vector<F> &_p);


	private:
		void vector(T _x);
		template<class... Ts>
		void vector(T _x, Ts... _xs);
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
	template<class ...Ts>
	inline Vector<T>::Vector(Ts ..._xs) : DOX(sizeof...(Ts)) {
		vector(_xs...);
	}


	template<class T>
	inline void Vector<T>::vector(T _x) {
		this->x.push_back(_x);
	}


	template<class T>
	template<class ...Ts>
	inline void Vector<T>::vector(T _x, Ts ..._xs) {
		this->x.push_back(_x);
		vector(_xs...);
	}


	template<class F>
	std::ostream & operator<<(std::ostream & _out, const Vector<F>& _p) {
		for (auto xi : _p.x) {
			_out << xi << "\t";
		}
		_out << std::endl;
		return _out;
	}
}