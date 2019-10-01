//*****************************************************************************
//Title		:PANSFEM2/LinearAlgebra/Models/Point.h
//Author	:Tanabe Yuta
//Date		:2019/10/01
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace PANSFEM2 {
	template<class T>
	class Point
	{
	public:
		Point();
		~Point();
		Point(std::vector<T> _x);


		const int DOX;


		template<class F>
		friend std::ostream& operator << (std::ostream &_out, const Point<F> &_p);	//streamÇ…èoóÕ


		std::vector<T> x;					//ç¿ïWíl
	};


	template<class T>
	inline Point<T>::Point() : DOX(0) {}


	template<class T>
	inline Point<T>::~Point() {}


	template<class T>
	inline Point<T>::Point(std::vector<T> _x) : DOX(_x.size()) {
		x = _x;
	}


	template<class F>
	inline std::ostream & operator<<(std::ostream & _out, const Point<F>& _p) {
		for (auto xi : _p.x) {
			_out << xi << "\t";
		}
		_out << std::endl;
		return _out;
	}
}