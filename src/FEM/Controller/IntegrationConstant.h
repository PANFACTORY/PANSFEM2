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
	const std::vector<std::vector<T> > GP3D8 = { 
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
	const std::vector<std::vector<T> > GW3D8 = { 
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