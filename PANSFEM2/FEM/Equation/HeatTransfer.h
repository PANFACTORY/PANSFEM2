//*****************************************************************************
//Title		:PANSFEM2/FEM/Equation/HeatTransfer.h
//Author	:Tanabe Yuta
//Date		:2019/10/03
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Vector.h"
#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
	//******************************�񎟌��M�`���̔M�`���}�g���N�X�𐶐�******************************
	template<class T>
	void HeatTransferTri(LILCSR<T>& _K, std::vector<Vector<T> >& _nodes, std::vector<int>& _nodestoelement, T _alpha) {

	}


	//******************************�񎟌��M�`���̔M�e�ʃ}�g���N�X�𐶐�******************************
	template<class T>
	void HeatCapacityTri(LILCSR<T>& _C, std::vector<Vector<T> >& _nodes, std::vector<int>& _nodestoelement, T _alpha) {

	}
}