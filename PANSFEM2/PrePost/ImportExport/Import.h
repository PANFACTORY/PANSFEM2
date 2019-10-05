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

		//.....一行読み飛ばす.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....一行分読み込む.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....節点IDを確認.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....節点の値を読み込む.....
				std::vector<double> xs;
				while(std::getline(sbuf, str, ',')) {
					xs.push_back(stod(str));
				}

				//.....節点を追加.....
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

		//.....一行読み飛ばす.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....一行分読み込む.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....要素IDを確認.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....要素内の節点IDを読み込む.....
				std::vector<int> element;
				while (std::getline(sbuf, str, ',')) {
					element.push_back(stoi(str));
				}

				//.....節点を追加.....
				_elements.push_back(element);
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportSystemIDFromCSV********************

}