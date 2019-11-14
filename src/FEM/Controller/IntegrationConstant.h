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
	//********************GaussIntegrationConstant2D1Points********************
	template<class T>
	class Gauss1Triangle{
public:
		static const int N = 1;								//Number of integration point
		static const std::vector<std::vector<T> > Points;	//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<std::vector<T> > Gauss1Triangle<T>::Points = { 
													{ 1.0 / 3.0, 1.0 / 3.0 },
												};


	template<class T>
	const std::vector<std::vector<T> > Gauss1Triangle<T>::Weights = { 
													{ 1.0 / sqrt(2.0), 1.0 / sqrt(2.0) },
												};
	

	//********************GaussIntegrationConstant2D4Points********************
	template<class T>
	class Gauss4Square{
public:
		static const int N = 4;								//Number of integration point
		static const std::vector<std::vector<T> > Points;	//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<std::vector<T> > Gauss4Square<T>::Points = { 
													{ -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0), -1.0 / sqrt(3.0) },
													{ -1.0 / sqrt(3.0),  1.0 / sqrt(3.0) },
													{  1.0 / sqrt(3.0),  1.0 / sqrt(3.0) }
												};


	template<class T>
	const std::vector<std::vector<T> > Gauss4Square<T>::Weights = { 
													{ 1.0, 1.0 },
													{ 1.0, 1.0 },
													{ 1.0, 1.0 },
													{ 1.0, 1.0 }
												};
	

	//********************GaussIntegrationConstant2D9Points********************
	template<class T>
	class Gauss9Square{
public:
		static const int N = 9;								//Number of integration point
		static const std::vector<std::vector<T> > Points;	//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<std::vector<T> > Gauss9Square<T>::Points = { 
													{ -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
													{ 0.0,				-sqrt(3.0 / 5.0) },
													{  sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
													{ -sqrt(3.0 / 5.0), 0.0				 },
													{ 0.0, 				0.0				 },
													{  sqrt(3.0 / 5.0), 0.0				 },
													{ -sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0) },
													{ 0.0,				 sqrt(3.0 / 5.0) },
													{  sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0) }
												};


	template<class T>
	const std::vector<std::vector<T> > Gauss9Square<T>::Weights = { 
													{ 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0 },
													{ 8.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0 }
												};


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


	//********************GaussIntegrationConstant27D8Points********************
	template<class T>
	class Gauss27Cubic{
public:
		static const int N = 27;							//Number of integration point
		static const std::vector<std::vector<T> > Points;	//Cordinate values of integration point
		static const std::vector<std::vector<T> > Weights;	//Weights of integration point
	};


	template<class T>
	const std::vector<std::vector<T> > Gauss27Cubic<T>::Points = {
													{ -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
													{  0.0,             -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
													{ sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
													{ -sqrt(3.0 / 5.0), 0.0,              -sqrt(3.0 / 5.0) },
													{ 0.0,              0.0,              -sqrt(3.0 / 5.0) },
													{ sqrt(3.0 / 5.0),  0.0,              -sqrt(3.0 / 5.0) },
													{ -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0) },
													{ 0.0,              sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0) },
													{ sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0) },
													{ -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0.0              },
													{ 0.0,              -sqrt(3.0 / 5.0), 0.0              },
													{ sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0), 0.0              },
													{ -sqrt(3.0 / 5.0), 0.0,              0.0              },
													{ 0.0,              0.0,              0.0              },
													{ sqrt(3.0 / 5.0),  0.0,              0.0              },
													{ -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0),  0.0              },
                                                    { 0.0,              sqrt(3.0 / 5.0),  0.0              },
													{ sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0),  0.0              },
													{ -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)  },
													{ 0.0,              -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)  },
													{ sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)  },
													{ -sqrt(3.0 / 5.0), 0.0,              sqrt(3.0 / 5.0)  },
													{ 0.0,              0.0,              sqrt(3.0 / 5.0)  },
													{ sqrt(3.0 / 5.0),  0.0,              sqrt(3.0 / 5.0)  },
													{ -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0)  },
													{ 0.0,              sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0)  },
													{ sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0),  sqrt(3.0 / 5.0)  }
												};


	template<class T>
	const std::vector<std::vector<T> > Gauss27Cubic<T>::Weights = { 
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
                                                    { 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0 },
													{ 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
                                                    { 8.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
                                                    { 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 },
													{ 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 }
												};
}