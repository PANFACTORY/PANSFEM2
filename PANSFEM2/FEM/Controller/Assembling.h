//*****************************************************************************
//Title		:PANSFEM2/FEM/Controller/Assembling.h
//Author	:Tanabe Yuta
//Date		:2019/10/04
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/LILCSR.h"


namespace PANSFEM2 {
	//********************Assembling K Matrix from Ke Matrix********************
	template<class T>
	void Assembling(LILCSR<T>& _K, std::vector<std::vector<T> >& _Ke, std::vector<int>& _element, std::vector<std::pair<int, int> >& _systemindex) {
		//.....アセンブリング.....
		int ei = 0;
		for (auto ni : _element) {
			for (int si = _systemindex[ni].first; si < _systemindex[ni].first + _systemindex[ni].second; si++) {
				int ej = 0;
				for (auto nj : _element) {
					for (int sj = _systemindex[nj].first; sj < _systemindex[nj].first + _systemindex[nj].second; sj++) {
						_K.set(si, sj, _K.get(si, sj) + _Ke[ei][ej]);
						ej++;
					}
				}
				ei++;
			}
		}
	}
}