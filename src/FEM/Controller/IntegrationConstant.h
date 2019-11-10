//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/IntegrationConstant.h
//Author	:Tanabe Yuta
//Date		:2019/10/12
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <vector>


namespace PANSFEM2 {
	//********************GaussIntegrationConstant3D8Points********************
	template<class T>
	class Gauss8Cubic{
public:
		static const int N = 8;								//Number of integration point
		static const std::vector<std::vector<T> > Points;	//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<std::vector<T> > Gauss8Cubic<T>::Points = { 
													{ -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0),  1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{ -1.0 / sqrt(3.0),  1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{ -1.0 / sqrt(3.0), -1.0 / sqrt(3.0),  1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0), -1.0 / sqrt(3.0),  1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0),  1.0 / sqrt(3.0),  1.0 / sqrt(3.0) },
													{ -1.0 / sqrt(3.0),  1.0 / sqrt(3.0),  1.0 / sqrt(3.0) } 
												};


	template<class T>
	const std::vector<std::vector<T> > Gauss8Cubic<T>::Weights = { 
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 },
													{ 1.0, 1.0, 1.0 } 
												};
}