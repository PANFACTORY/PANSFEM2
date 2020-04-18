//*****************************************************************************
//Title		:PANSFEM2/PrePost/Import/ImportFromCSV.h
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

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of a node.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....Read values of a node.....
				std::vector<T> x;
				while(std::getline(sbuf, str, ',')) {
					x.push_back(stod(str));
				}

				//.....Add a node.....
				_nodes.push_back(Vector<T>(x));
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

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of an element.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				//.....Get ids of nodes in an element.....
				std::vector<int> element;
				while (std::getline(sbuf, str, ',')) {
					element.push_back(stoi(str));
				}

				//.....Add an element.....
				_elements.push_back(element);
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportDirichletConditionFromCSV********************
	template<class T>
	bool ImportDirichletFromCSV(std::vector<std::pair<std::pair<int, int>, T> >& _ufixed, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Dirichlet Condition file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of a node.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				int idn = stoi(str);
				for(int i = 0; std::getline(sbuf, str, ','); i++) {
					if(str != "free") {
						_ufixed.push_back(std::make_pair(std::make_pair(idn, i), stod(str)));
					}
				}
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportNeumannConditionFromCSV********************
	template<class T>
	bool ImportNeumannFromCSV(std::vector<std::pair<std::pair<int, int>, T> >& _qfixed, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Dirichlet Condition file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of a node.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				int idn = stoi(str);
				for(int i = 0; std::getline(sbuf, str, ','); i++) {
					if(str != "free") {
						_qfixed.push_back(std::make_pair(std::make_pair(idn, i), stod(str)));
					}
				}
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportInitialConditionFromCSV********************
	template<class T>
	bool ImportInitialFromCSV(std::vector<Vector<T> >& _u, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Initial Condition file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of a node.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				int idn = stoi(str);
				for (int i = 0; i < _u[idn].SIZE(); i++) {
					std::getline(sbuf, str, ',');
					if (str != "free") {
						_u[idn](i) = stod(str);
					}
				}
			}
		}

		ifs.close();

		return true;
	}


	//********************ImportPeriodicBoundaryConditionFromCSV*******************
	bool ImportPeriodicFromCSV(std::vector<std::pair<int, int> >& _ufixed, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Periodic Boundary Condition file " << _fname << " open error!" << std::endl;
			return false;
		}

		//.....Pass a line.....
		std::string str0;
		std::getline(ifs, str0);

		while (!ifs.eof()) {
			//.....Read a line.....
			std::string buf;
			ifs >> buf;
			std::istringstream sbuf(buf);
			std::string str;

			//.....Get id of a node.....
			std::getline(sbuf, str, ',');
			if (!str.empty()) {
				int masterid = stoi(str);
				std::getline(sbuf, str, ',');
				int slaveid = stoi(str);
				_ufixed.push_back(std::make_pair(masterid, slaveid));
			}
		}

		ifs.close();

		return true;
	}
}