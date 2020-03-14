//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/NewtonCotesIntegration.h
//Author	:Tanabe Yuta
//Date		:2020/03/14
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <vector>
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************NewtonCotesIntegrationConstant1Line********************
	template<class T>
	class NewtonCotes1Line{
public:
		static const int N = 1;								//Number of integration point
		static const std::vector<Vector<T> > Points;		//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<Vector<T> > NewtonCotes1Line<T>::Points = { 
													{ 0.0 },
												};


	template<class T>
	const std::vector<std::vector<T> > NewtonCotes1Line<T>::Weights = { 
													{ 2.0 },
												};


    //********************NewtonCotesIntegrationConstant3Line********************
	template<class T>
	class NewtonCotes3Line{
public:
		static const int N = 3;								//Number of integration point
		static const std::vector<Vector<T> > Points;		//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<Vector<T> > NewtonCotes3Line<T>::Points = { 
													{ -1.0 },
                                                    { 0.0 },
                                                    { 1.0 },
												};


	template<class T>
	const std::vector<std::vector<T> > NewtonCotes3Line<T>::Weights = { 
													{ 2.0/6.0 },
                                                    { 8.0/6.0 },
                                                    { 2.0/6.0 },
												};


    //********************NewtonCotesIntegrationConstant5Line********************
	template<class T>
	class NewtonCotes5Line{
public:
		static const int N = 5;								//Number of integration point
		static const std::vector<Vector<T> > Points;		//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<Vector<T> > NewtonCotes5Line<T>::Points = { 
													{ -1.0 },
                                                    { -0.5 },
                                                    { 0.0 },
                                                    { 0.5 },
                                                    { 1.0 },
												};


	template<class T>
	const std::vector<std::vector<T> > NewtonCotes5Line<T>::Weights = { 
													{ 14.0/90.0 },
                                                    { 64.0/90.0 },
                                                    { 24.0/90.0 },
                                                    { 64.0/90.0 },
                                                    { 14.0/90.0 },
												};


    //********************NewtonCotesIntegrationConstant7Line********************
	template<class T>
	class NewtonCotes7Line{
public:
		static const int N = 7;								//Number of integration point
		static const std::vector<Vector<T> > Points;		//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<Vector<T> > NewtonCotes7Line<T>::Points = { 
													{ -1.0 },
                                                    { -2.0/3.0 },
                                                    { -1.0/3.0 },
                                                    { 0.0 },
                                                    { 1.0/3.0 },
                                                    { 2.0/3.0 },
                                                    { 1.0 },
												};


	template<class T>
	const std::vector<std::vector<T> > NewtonCotes7Line<T>::Weights = { 
													{ 82.0/840.0 },
                                                    { 432.0/840.0 },
                                                    { 54.0/840.0 },
                                                    { 544.0/840.0 },
                                                    { 54.0/840.0 },
                                                    { 432.0/840.0 },
                                                    { 82.0/840.0 },
												};
}