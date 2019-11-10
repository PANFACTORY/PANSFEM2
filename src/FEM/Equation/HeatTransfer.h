//*****************************************************************************
//Title		:src/FEM/Equation/HeatTransfer.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//******************************Heat transfer matrix******************************
	template<class T>
	HeatTransferTri(Matrix<T>& _Ke, std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _alpha, T _t) {
		//.....Get element matrix.....
		_Ke = Matrix<T>(_element.size(), _element.size());

		//.....Get space of element.....
		T A = 0.5*((_nodes[_element[0]](0) - _nodes[_element[2]](0))*(_nodes[_element[1]](1) - _nodes[_element[2]](1)) - (_nodes[_element[2]](1) - _nodes[_element[0]](1))*(_nodes[_element[2]](0) - _nodes[_element[1]](0)));

		//.....Make B matrix.....
		Matrix<T> B = Matrix<T>(2, 3);
		B(0, 0) = _nodes[_element[1]](1) - _nodes[_element[2]](1);	B(0, 1) = _nodes[_element[2]](1) - _nodes[_element[0]](1);	B(0, 2) = _nodes[_element[0]](1) - _nodes[_element[1]](1);
		B(1, 0) = _nodes[_element[2]](0) - _nodes[_element[1]](0);	B(1, 1) = _nodes[_element[0]](0) - _nodes[_element[2]](0);	B(1, 2) = _nodes[_element[1]](0) - _nodes[_element[0]](0);
		B /= (2.0*A);

		_Ke += _alpha * B.Transpose()*B*_t;
	}


	//******************************Heat capacity matrix******************************
	template<class T>
	HeatCapacityTri(Matrix<T>& _Ce, std::vector<Vector<T> >& _nodes, std::vector<int>& _element, T _rho, T _c, T _t) {
		//.....Get element matrix.....
		_Ce = Matrix<T>(_element.size(), _element.size());
		
		//.....Get space of element.....
		T A = 0.5*((_nodes[_element[0]][0] - _nodes[_element[2]][0])*(_nodes[_element[1]][1] - _nodes[_element[2]][1]) - (_nodes[_element[2]][1] - _nodes[_element[0]][1])*(_nodes[_element[2]][0] - _nodes[_element[1]][0]));

		//.....Make C matrix.....
		_Ce(0, 0) = 1.0 / 6.0;	_Ce(0, 1) = 1.0 / 12.0;	_Ce(0, 2) = 1.0 / 12.0;
		_Ce(1, 0) = 1.0 / 12.0;	_Ce(1, 1) = 1.0 / 6.0;	_Ce(1, 2) = 1.0 / 12.0;
		_Ce(2, 0) = 1.0 / 12.0;	_Ce(2, 1) = 1.0 / 12.0;	_Ce(2, 2) = 1.0 / 6.0;
		_Ce *= (_rho * _c*A*_t);
	}
}