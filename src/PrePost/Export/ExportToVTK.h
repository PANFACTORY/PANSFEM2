//*****************************************************************************
//Title		:PANSFEM2/PrePost/Export/ExportToVTK.h
//Author	:Tanabe Yuta
//Date		:2019/10/06
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <fstream>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************Make Headder********************
	void MakeHeadderToVTK(std::ofstream& _fout) {
		_fout << "# vtk DataFile Version 4.1\n";
		_fout << "vtk output\n";
		_fout << "ASCII\n";
		_fout << "DATASET UNSTRUCTURED_GRID\n";
	}


	//********************Add Points********************
	template<class T>
	void AddPointsToVTK(std::vector<Vector<T> > _nodes, std::ofstream& _fout) {
		_fout << "\nPOINTS\t" << _nodes.size() << "\tfloat\n";
		for (auto node : _nodes) {
			int i = 0;
			for(; i < node.SIZE(); i++){
				_fout << node(i) << "\t";
			}
			for(; i < 3; i++){
				_fout << T() << "\t";
			}
			_fout << std::endl;
		}
	}


	//********************Add Elements*******************
	void AddElementToVTK(std::vector<std::vector<int> > _elements, std::ofstream& _fout) {
		int datanum = 0;
		for (auto element : _elements) {
			datanum += element.size() + 1;
		}

		_fout << "\nCELLS " << _elements.size() << "\t" << datanum << "\n";
		for (auto element : _elements) {
			_fout << element.size() << "\t";
			for (auto node : element) {
				_fout << node << "\t";
			}
			_fout << std::endl;
		}
	}


	//********************Add Element Types********************
	void AddElementTypes(std::vector<int> _elementtypes, std::ofstream& _fout) {
		_fout << "\nCELL_TYPES\t" << _elementtypes.size() << "\n";
		for (auto elementtype : _elementtypes) {
			_fout << elementtype << "\n";
		}
	}


	//********************Add Point Scalers********************
	template<class T>
	void AddPointScalers(std::vector<T> _values, std::string _symbol, std::ofstream& _fout) {
		_fout << "\nPOINT_DATA\t" << _values.size() << "\n";
		_fout << "SCALARS " << _symbol << " float\n";
		_fout << "LOOKUP_TABLE default\n";

		for (auto value : _values) {
			_fout << value << std::endl;
		}
	}


	//********************Add Point Vectors********************
	template<class T>
	void AddPointVectors(std::vector<Vector<T> > _values, std::string _symbol, std::ofstream& _fout) {
		_fout << "\nPOINT_DATA\t" << _values.size() << "\n";
		_fout << "VECTORS " << _symbol << " float\n";

		for (auto value : _values) {
			int i = 0;
			for(; i < value.SIZE(); i++){
				_fout << value(i) << "\t";
			}
			for(; i < 3; i++){
				_fout << T() << "\t";
			}
			_fout << std::endl;
		}
	}


	//********************Add Point Tensors********************


	//********************Add Element Scalers********************
	template<class T>
	void AddElementScalers(std::vector<T> _values, std::string _symbol, std::ofstream& _fout) {
		_fout << "\nCELL_DATA\t" << _values.size() << "\n";
		_fout << "SCALARS " << _symbol << " float\n";
		_fout << "LOOKUP_TABLE default\n";

		for (auto value : _values) {
			_fout << value << std::endl;
		}
	}


	//********************Add Element Vectors********************
	//********************Add Element Tensors********************
}