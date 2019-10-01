//*****************************************************************************
//Title		:PANSFEM2/FEM/Models/Mesh.h
//Author	:Tanabe Yuta
//Date		:2019/10/01
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>
#include <vector>


#include "../../LinearAlgebra/Models/Point.h"


namespace PANSFEM2 {
	template<class T>
	class Mesh
	{
	public:
		Mesh();
		~Mesh();
		Mesh(int _DOX);


		const int DOX;		//�Ɨ��ϐ�����


		void setpnode(Point<T>* _pnode);					//�ߓ_���w���|�C���^��ǉ�
		void setelement(std::vector<int> _nodestoelement);	//�ߓ_�\�v�f�֌W��ǉ�


	private:
		std::vector<Point<T>*> pnodes;						//�ߓ_���w���|�C���^
		std::vector<std::vector<int> > nodestoelements;		//�ߓ_�\�v�f�֌W
	};


	template<class T>
	inline Mesh<T>::Mesh() : DOX(0) {}


	template<class T>
	inline Mesh<T>::~Mesh() {}


	template<class T>
	inline Mesh<T>::Mesh(int _DOX) : DOX(_DOX) {}


	template<class T>
	inline void Mesh<T>::setpnode(Point<T>* _pnode)	{
		assert(_pnode->DOX == this->DOX);
		this->pnodes.push_back(_pnode);
	}


	template<class T>
	inline void Mesh<T>::setelement(std::vector<int> _nodestoelement) {
		this->nodestoelements.push_back(_nodestoelement);
	}
}