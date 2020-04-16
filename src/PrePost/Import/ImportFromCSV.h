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


	//********************ImportSystemIDFromCSV********************
	bool ImportFieldFromCSV(std::vector<int>& _field, int& _currentindex, int _nodessize, std::string _fname) {
		_field = std::vector<int>(_nodessize, 0);

		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Field file " << _fname << " open error!" << std::endl;
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
				_field[idn] = _currentindex;
				
				std::getline(sbuf, str, ',');
				int size = stoi(str);
				_currentindex += size;
			}
		}

		_field.push_back(_currentindex);

		ifs.close();

		return true;
	}


	//********************ImportDirichletConditionFromCSV********************
	template<class T>
	bool ImportDirichletFromCSV(std::vector<int>& _isufixed, std::vector<T>& _ufixed, std::vector<int> _field, std::string _fname) {
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
				for (int i = _field[idn]; i < _field[idn + 1]; i++) {
					std::getline(sbuf, str, ',');
					if (str != "free") {
						_isufixed.push_back(i);
						_ufixed.push_back(stod(str));
					}
				}
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
	bool ImportNeumannFromCSV(std::vector<int>& _isqfixed, std::vector<T>& _qfixed, std::vector<int> _field, std::string _fname) {
		std::ifstream ifs(_fname);

		if (!ifs.is_open()) {
			std::cout << "Neumann Condition file " << _fname << " open error!" << std::endl;
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
				for (int i = _field[idn]; i < _field[idn + 1]; i++) {
					std::getline(sbuf, str, ',');
					if (str != "free") {
						_isqfixed.push_back(i);
						_qfixed.push_back(stod(str));
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
	bool ImportPeriodicFromCSV(std::vector<int>& _ismasterfixed, std::vector<int>& _isslavefixed, std::vector<int> _field, std::string _fname) {
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

				//.....Convert id to fix date.....
				if((_field[masterid + 1] - _field[masterid]) != (_field[slaveid + 1] - _field[slaveid])) {
					return false;
				}
				for(int i = 0; _field[masterid] + i < _field[masterid + 1]; i++) {
					_ismasterfixed.push_back(_field[masterid] + i);
					_isslavefixed.push_back(_field[slaveid] + i);
				}
			}
		}

		ifs.close();

		return true;
	}
}