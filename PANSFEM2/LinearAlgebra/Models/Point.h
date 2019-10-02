//*****************************************************************************
//Title		:PANSFEM2/LinearAlgebra/Models/Point.h
//Author	:Tanabe Yuta
//Date		:2019/10/02
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <vector>


namespace PANSFEM2 {
	template<class T>
	class Point
	{
	public:
		Point();
		~Point();
		template<class... Ts>
		Point(Ts... _xs);


		const int DOX;
		std::vector<T> x;


		template<class F>
		friend std::ostream& operator << (std::ostream &_out, const Point<F> &_p);


	private:
		void point(T _x);
		template<class... Ts>
		void point(T _x, Ts... _xs);
	};


	template<class T>
	inline Point<T>::Point() : DOX(0) {}


	template<class T>
	inline Point<T>::~Point() {}

	
	template<class T>
	template<class ...Ts>
	inline Point<T>::Point(Ts ..._xs) : DOX(sizeof...(Ts)) {
		point(_xs...);
	}


	template<class T>
	inline void Point<T>::point(T _x) {
		this->x.push_back(_x);
	}


	template<class T>
	template<class ...Ts>
	inline void Point<T>::point(T _x, Ts ..._xs) {
		this->x.push_back(_x);
		point(_xs...);
	}


	template<class F>
	std::ostream & operator<<(std::ostream & _out, const Point<F>& _p) {
		for (auto xi : _p.x) {
			_out << xi << "\t";
		}
		_out << std::endl;
		return _out;
	}
}