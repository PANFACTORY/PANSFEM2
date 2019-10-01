//*****************************************************************************
//Title		:PANSFEM2/FEM/Models/Field.h
//Author	:Tanabe Yuta
//Date		:2019/10/01
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "Mesh.h"


namespace PANSFEM2 {
	template<class T>
	class Field
	{
	public:
		Field();
		~Field();
		Field(Mesh<T>* _pmesh, int _DOU, int _DOS);


		const int DOU;		//従属変数次元
		const int DOS;		//従属変数ステップ数


		Mesh<T>* pmesh;							//Meshを指すポインタ
		std::vector <std::vector<T> > values;	//従属変数値
	};


	template<class T>
	inline Field<T>::Field() : DOU(0), DOS(0) {}


	template<class T>
	inline Field<T>::~Field() {}


	template<class T>
	inline Field<T>::Field(Mesh<T>* _pmesh, int _DOU, int _DOS) : DOU(_DOU), DOS(_DOS) {
		this->pmesh = _pmesh;
	}
}