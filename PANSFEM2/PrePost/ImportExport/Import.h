//*****************************************************************************
//Title		:PANSFEM2/PrePost/ImportExport/Import.h
//Author	:Tanabe Yuta
//Date		:2019/10/01
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
	//********************ImportNodesFromCSV********************
	template<class T>
	bool ImportNodesFromCSV(std::vector<Vector<T> >& _nodes, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Node file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....��s�ǂݔ�΂�.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....��s���ǂݍ���.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....�ߓ_ID���m�F.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....�ߓ_�̒l��ǂݍ���.....
				std::vector<double> xs;
				while(std::getline(sbuf, str, ',')) {
					xs.push_back(stod(str));
				}

				//.....�ߓ_��ǉ�.....
				_nodes.push_back(Vector<T>(xs));
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportElementFromCSV********************
	bool ImportElementsFromCSV(std::vector<std::vector<int> >& _elements, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Element file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....��s�ǂݔ�΂�.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....��s���ǂݍ���.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....�v�fID���m�F.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....�v�f���̐ߓ_ID��ǂݍ���.....
				std::vector<int> element;
				while (std::getline(sbuf, str, ',')) {
					element.push_back(stoi(str));
				}

				//.....�ߓ_��ǉ�.....
				_elements.push_back(element);
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportSystemIDFromCSV********************

}