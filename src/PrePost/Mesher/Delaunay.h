//*****************************************************************************
//  Title		:   src/PrePost/Mesher/Delaunay.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/04/27
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <array>


#include "../../LinearAlgebra/Models/Vector.h"


#define SELFEPS0 1.0e-10
#define SELFEPS1 1.0e-10


namespace PANSFEM2{
    //********************Delaunay Boundary class********************
    class Boundary {
public:
        Boundary(std::vector<int>& _nodelists, bool _type);
		~Boundary();
		
		std::vector<int> nodelists;		//list of nodes on boundary
		bool type;						//type of boundary 

		int order(int _nodenum);		//get order id on boundary
	};


	Boundary::Boundary(std::vector<int>& _nodelists, bool _type){
		this->nodelists = _nodelists;
		this->type = _type;
	}


	Boundary::~Boundary(){}


	int Boundary::order(int _nodenum) {
		for (int i = 0; i < this->nodelists.size(); i++) {
			if (this->nodelists[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}


    //********************Delaunay Node class********************
    template<class T>
	class Node {
public:
		Node();
		~Node();
        Node(Vector<T> _x);
		Node(T _x, T _y);

		Vector<T> x;			//cordinate vallue of Node
		bool isset;				//is set						
		bool isonboundary;		//is on boundary

		T distance(Node<T> _node);							//get distance between other node
		T vecpro(Node<T> _node0, Node<T> _node1);		//get innerproduct
		T innpro(Node<T> _node0, Node<T> _node1);		//get vectorproduct
	};


	template<class T>
	Node<T>::Node(){
		this->x = { T(), T() };
		this->isset = false;
		this->isonboundary = false;
	}


	template<class T>
	Node<T>::~Node(){}


    template<class T>
    Node<T>::Node(Vector<T> _x) {
        this->x = _x;
        this->isset = false;
		this->isonboundary = false;
    }


	template<class T>
	Node<T>::Node(T _x, T _y){
		this->x = { _x, _y };
		this->isset = false;
		this->isonboundary = false;
	}


	template<class T>
	T Node<T>::distance(Node<T> _node) {
		return (this->x - _node.x).Norm();
	}


	template<class T>
	T Node<T>::vecpro(Node<T> _node0, Node<T> _node1) {
		return (_node0.x(0) - this->x(0))*(_node1.x(1) - this->x(1)) - (_node0.x(1) - this->x(1))*(_node1.x(0) - this->x(0));
	}


	template<class T>
	T Node<T>::innpro(Node<T> _node0, Node<T> _node1) {
		return (_node0.x - this->x)*(_node1.x - this->x);
	}


    //********************Delaunay Element class********************
    template<class T>
	class Element {
public:
		Element();
		~Element();
		Element(int _node0, int _node1, int _node2);

		std::array<int, 3> nodes;								//id of nodes
		std::array<int, 3> neighbors;							//id of neighbor element
		std::array<T, 3> angles;								//angle of each corner
		std::array<bool, 3> sides;								//is edge on boundary
		bool active;											//is in boundary
		bool check;												//has already checked

		void setnode(int _node0, int _node1, int _node2);					//set node id
		void setneighbor(int _neighbor0, int _neighbor1, int _neighbor2);	//set neighbor id
		void setside(bool _side0, bool _side1, bool _side2);				//set edge is or not on boundary
		void copy(const Element<T>& _originalelement);						//copy element
		void getangle(std::vector<Node<T> >& _node);						//get angle

		int inouton(int _nodenum, std::vector<Node<T> >& _nodes);			//get location of node
		int oppositenode(int _elementname);									//get id of opposite node 

		T space(std::vector<Node<T> >& _node);								//get space of element
		int nodeorder(int _nodenum);										//get node id in element 
private: 
	};


	template<class T>
	Element<T>::Element(){
		this->nodes[0] = -1;
		this->nodes[1] = -1;
		this->nodes[2] = -1;

		this->sides[0] = false;
		this->sides[1] = false;
		this->sides[2] = false;

		this->neighbors[0] = -1;
		this->neighbors[1] = -1;
		this->neighbors[2] = -1;

		this->active = true;
		this->check = false;
	}


	template<class T>
	Element<T>::~Element(){}


	template<class T>
	Element<T>::Element(int _node0, int _node1, int _node2){
		this->nodes[0] = _node0;
		this->nodes[1] = _node1;
		this->nodes[2] = _node2;

		this->sides[0] = false;				
		this->sides[1] = false;
		this->sides[2] = false;

		this->neighbors[0] = -1;				
		this->neighbors[1] = -1;
		this->neighbors[2] = -1;

		this->active = true;
		this->check = false;
	}


	template<class T>
	void Element<T>::setnode(int _node0, int _node1, int _node2) {
		this->nodes[0] = _node0;
		this->nodes[1] = _node1;
		this->nodes[2] = _node2;
	}


	template<class T>
	void Element<T>::setneighbor(int _neighbor0, int _neighbor1, int _neighbor2) {
		this->neighbors[0] = _neighbor0;
		this->neighbors[1] = _neighbor1;
		this->neighbors[2] = _neighbor2;
	}


	template<class T>
	void Element<T>::setside(bool _side0, bool _side1, bool _side2) {
		this->sides[0] = _side0;
		this->sides[1] = _side1;
		this->sides[2] = _side2;
	}


	template<class T>
	void Element<T>::copy(const Element<T>& _originalelement) {
		this->nodes = _originalelement.nodes;
		this->neighbors = _originalelement.neighbors;
		this->sides = _originalelement.sides;
		this->angles = _originalelement.angles;
	}


	//*****************************************************************************
	//	return -(i+1)	�Fi�Ԗڂ̕ӂ̊O��
	//	return i+1		�Fi�Ԗڂ̕ӏ�
	//	return 0		�F�O�p�`����
	//*****************************************************************************
	template<class T>
	int Element<T>::inouton(int _nodenum, std::vector<Node<T> >& _nodes) {
		T vecpro0 = _nodes[this->nodes[0]].vecpro(_nodes[this->nodes[1]], _nodes[_nodenum]);
		T vecpro1 = _nodes[this->nodes[1]].vecpro(_nodes[this->nodes[2]], _nodes[_nodenum]);
		T vecpro2 = _nodes[this->nodes[2]].vecpro(_nodes[this->nodes[0]], _nodes[_nodenum]);

		T vecpro3 = _nodes[this->nodes[0]].vecpro(_nodes[this->nodes[2]], _nodes[_nodenum]);
		T vecpro4 = _nodes[this->nodes[1]].vecpro(_nodes[this->nodes[0]], _nodes[_nodenum]);
		T vecpro5 = _nodes[this->nodes[2]].vecpro(_nodes[this->nodes[1]], _nodes[_nodenum]);

		if (vecpro2 > T() && (vecpro0 < T() || (fabs(vecpro0) <= SELFEPS0 && vecpro5 > T()))) {
			return -3;
		}
		
		else if (vecpro0 > T() && (vecpro1 < T() || (fabs(vecpro1) <= SELFEPS0 && vecpro3 > T()))) {
			return -1;
		}
		
		else if (vecpro1 > T() && (vecpro2 < T() || (fabs(vecpro2) <= SELFEPS0 && vecpro4 > T()))) {
			return -2;
		}
		
		else if (fabs(vecpro0) <= SELFEPS1) {
			return 3;
		}
		
		else if (fabs(vecpro1) <= SELFEPS1) {
			return 1;
		}
		
		else if (fabs(vecpro2) <= SELFEPS1) {
			return 2;
		}
		return 0;
	}


	template<class T>
	int Element<T>::nodeorder(int _nodenum) {
		for (int i = 0; i < 3; i++) {
			if (this->nodes[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	int Element<T>::oppositenode(int _elementname) {
		for (int i = 0; i < 3; i++) {
			if (this->neighbors[i] == _elementname) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	void Element<T>::getangle(std::vector<Node<T> >& _nodes) {
		for (int i = 0; i < 3; i++) {
			this->angles[i] = 180.0*acos(_nodes[this->nodes[i]].innpro(_nodes[this->nodes[(i + 1) % 3]], _nodes[this->nodes[(i + 2) % 3]]) / (_nodes[this->nodes[i]].distance(_nodes[this->nodes[(i + 1) % 3]]) * _nodes[this->nodes[i]].distance(_nodes[this->nodes[(i + 2) % 3]]))) / M_PI;
		}
	}


	template<class T>
	T Element<T>::space(std::vector<Node<T> >& _nodes) {
		return 0.5*_nodes[this->nodes[0]].vecpro(_nodes[this->nodes[1]], _nodes[this->nodes[2]]);
	}


	//********************Delaunay mesher class********************
	template<class T>
	class Delaunay {
public:
		Delaunay(std::vector<Vector<T> > _x, std::vector<std::vector<int> > _boundaries, T _maxsize);
		~Delaunay();

        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();
private:
		const int ADDITIONALNODENUM0 = 100000;

		std::vector<Node<T> > nodes;
		std::vector<Element<T> > elements;

		void getsupertriangle();
		void getboundary(Boundary _boundary);
		void deactivate(Boundary _boundary);
		void deletesupertriangle();
		void deleteelement();
		void sortelement();
		void getinternalelement(T _maxside);

		void getelementin(int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1);
		void getelementon(int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1);
		void swapping(std::vector<int>& _stack, int _nodenump1, int _nodenumm1);
	};


	template<class T>
	Delaunay<T>::Delaunay(std::vector<Vector<T> > _x, std::vector<std::vector<int> > _boundaries, T _maxsize) {
		//----------Generate nodes, elements and boundaries----------
		this->nodes = std::vector<Node<T> >();
        for(auto xi : _x) {
            this->nodes.push_back(Node<T>(xi));
        }

		this->elements = std::vector<Element<T> >();

        std::vector<Boundary> boundaries = std::vector<Boundary>();
        for(auto boundary : _boundaries) {
            boundaries.push_back(Boundary(boundary, true));
        }
		
		//----------Generate SuperTriangle----------
		this->getsupertriangle();

		//----------Generate Boundary----------
		for (auto boundary : boundaries) {
			this->getboundary(boundary);
		}
		for (auto boundary : boundaries) {
			this->deactivate(boundary);
		}

		//----------Delete needless Elements----------
		this->deletesupertriangle();
		this->deleteelement();

		//----------Sort Elements----------
		this->sortelement();

		//----------subdivide----------
		if (_maxsize > 0) {
			this->getinternalelement(_maxsize);
		}
	}


	template<class T>
	Delaunay<T>::~Delaunay() {}


	template<class T>
	void Delaunay<T>::getsupertriangle() {
		//Get distance maximam
		T rmax = T();
		Node<T> o;
		for (auto& node : this->nodes) {
			T tmpr = node.distance(o);
			if (rmax < tmpr) {
				rmax = tmpr;
			}
		}
		
		//Generate SuperTriangle
		for (int i = 0; i < 3; i++) {
			this->nodes.push_back(Node<T>(-2.0*rmax*sin(2.0*M_PI*i/3.0), 2.0*rmax*cos(2.0*M_PI*i/3.0)));
		}
		this->elements.push_back(Element<T>(this->nodes.size() - 3, this->nodes.size() - 2, this->nodes.size() - 1));
	}


	template<class T>
	void Delaunay<T>::getboundary(Boundary _boundary) {
		for (int i = 0; i < _boundary.nodelists.size(); i++) {
			//.....Add nodes on boundary into district.....
			if (this->nodes[_boundary.nodelists[i]].isset == false) {
				this->nodes[_boundary.nodelists[i]].isset = true;
				int nowtri = 0;
				if (this->elements.size() > 0) {
					nowtri = this->elements.size() - 1;
				}
				this->nodes[_boundary.nodelists[i]].isonboundary = true;

				//.....Search element in which node is.....
				for (int j = 0; j < this->elements.size(); j++) {
					int pos = this->elements[nowtri].inouton(_boundary.nodelists[i], this->nodes);
					//if not in or out
					if (pos < 0 || this->elements[nowtri].active == false) {
						if (this->elements[nowtri].neighbors[abs(pos) - 1] >= 0) {
							nowtri = this->elements[nowtri].neighbors[abs(pos) - 1];
						} else {
							std::cout << "Out of triangle Error!\n";
						}
					}
					//if in 
					else if (pos == 0) {
						if (i == 0) {
							getelementin(nowtri, _boundary.nodelists[i + 1], _boundary.nodelists[i], -2);
						} else if (i == _boundary.nodelists.size() - 1) {
							getelementin(nowtri, _boundary.nodelists[0], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						} else {
							getelementin(nowtri, _boundary.nodelists[i + 1], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						}
						break;
					}
					//if on
					else {
						if (i == 0) {
							getelementon(nowtri, pos - 1, _boundary.nodelists[i + 1], _boundary.nodelists[i], -2);
						} else if (i == _boundary.nodelists.size() - 1) {
							getelementon(nowtri, pos - 1, _boundary.nodelists[0], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						} else {
							getelementon(nowtri, pos - 1, _boundary.nodelists[i + 1], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						}
						break;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::deactivate(Boundary _boundary) {
		for (auto& element : this->elements) {
			if(element.check == false){
				//.....get order of node on boundary.....
				std::array<int, 3> nodeorders{ _boundary.order(element.nodes[0]), _boundary.order(element.nodes[1]), _boundary.order(element.nodes[2]) };
				
				//.....external boundary.....
				if (_boundary.type == true) {
					if (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0) {
						element.check = true;
						if ((nodeorders[0] < nodeorders[1] && nodeorders[1] < nodeorders[2])
							|| (nodeorders[1] < nodeorders[2] && nodeorders[2] < nodeorders[0])
							|| (nodeorders[2] < nodeorders[0] && nodeorders[0] < nodeorders[1])) {
							element.active = false;
						}
					}
				}
				//.....internal boundary.....
				else {
					//.....some nodes are on other boundary.....
					if ((nodeorders[0] < 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] < 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] < 0)
						|| (nodeorders[0] < 0 && nodeorders[1] < 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] < 0 && nodeorders[2] < 0)
						|| (nodeorders[0] < 0 && nodeorders[1] >= 0 && nodeorders[2] < 0)) {
						element.check = true;
					} else if (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0) {
						element.check = true;
						if ((nodeorders[0] < nodeorders[1] && nodeorders[1] < nodeorders[2])
							|| (nodeorders[1] < nodeorders[2] && nodeorders[2] < nodeorders[0])
							|| (nodeorders[2] < nodeorders[0] && nodeorders[0] < nodeorders[1])) {
							element.active = false;
						}
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::deletesupertriangle() {
		for (int i = this->elements.size() - 1; i >= 0; i--) {
			for (const auto& node : this->elements[i].nodes) {
				if (node == this->nodes.size() - 1 || node == this->nodes.size() - 2 || node == this->nodes.size() - 3) {
					this->elements[i].active = false;
					break;
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::deleteelement() {
		for (int i = this->elements.size() - 1; i >= 0; i--) {
			if (this->elements[i].active == false || this->elements[i].check == false) {
				this->elements.erase(this->elements.begin() + i);
			}
		}
	}


	template<class T>
	void Delaunay<T>::sortelement() {
		for (int i = 0; i < this->elements.size(); i++) {
			for (int j = 0; j < 3; j++) {
				this->elements[i].neighbors[j] = -1;
				for (int k = 0; k < this->elements.size(); k++) {
					if (i != k && ((this->elements[i].nodes[(j + 1)%3] == this->elements[k].nodes[0] && this->elements[i].nodes[(j + 2)%3] == this->elements[k].nodes[2])
						|| (this->elements[i].nodes[(j + 1)%3] == this->elements[k].nodes[1] && this->elements[i].nodes[(j + 2)%3] == this->elements[k].nodes[0])
						|| (this->elements[i].nodes[(j + 1)%3] == this->elements[k].nodes[2] && this->elements[i].nodes[(j + 2)%3] == this->elements[k].nodes[1]))) {
						this->elements[i].neighbors[j] = k;
						break;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::getelementin(int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		std::array<Element<T>, 3> tmpelement;

		for(int i = 0; i < 3; i++){
			tmpelement[i].sides[0] = this->elements[_nowtri].sides[i];
			if (this->elements[_nowtri].nodes[i] == _nodenumm1 || this->elements[_nowtri].nodes[i] == _nodenump1) {
				tmpelement[(i+1)%3].sides[1] = true;
				tmpelement[(i+2)%3].sides[2] = true;
			}
			tmpelement[i].setnode(_nodenum, this->elements[_nowtri].nodes[(i+1)%3], this->elements[_nowtri].nodes[(i+2)%3]);
			tmpelement[i].getangle(this->nodes);
		}
		tmpelement[0].setneighbor(this->elements[_nowtri].neighbors[0], this->elements.size(), this->elements.size() + 1);
		tmpelement[1].setneighbor(this->elements[_nowtri].neighbors[1], this->elements.size() + 1, _nowtri);
		tmpelement[2].setneighbor(this->elements[_nowtri].neighbors[2], _nowtri, this->elements.size());

		for (int k = 0; k < 2; k++) {
			int neighbor = this->elements[_nowtri].neighbors[1 + k];
			if (neighbor >= 0) {
				this->elements[neighbor].neighbors[this->elements[neighbor].oppositenode(_nowtri)] = this->elements.size() + k;
			}
		}

		stack.push_back(_nowtri);
		stack.push_back(this->elements.size());
		stack.push_back(this->elements.size() + 1);

		this->elements[_nowtri].copy(tmpelement[0]);
		this->elements.push_back(tmpelement[1]);
		this->elements.push_back(tmpelement[2]);

		this->swapping(stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void Delaunay<T>::getelementon(int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		int nownode = _pos;
		int neitri = this->elements[_nowtri].neighbors[nownode];
		
		//there is neighbor element
		if (neitri != -1 && this->elements[neitri].active == true) {
			int neinode = this->elements[neitri].oppositenode(_nowtri);
			std::array<Element<T>, 4> tmptri;			//0:nowtri	1:neitri

			tmptri[0].setnode(_nodenum, this->elements[_nowtri].nodes[nownode], this->elements[_nowtri].nodes[(nownode + 1)%3]);
			tmptri[0].setneighbor(this->elements[_nowtri].neighbors[(nownode + 2)%3], neitri, this->elements.size());
			tmptri[0].setside(this->elements[_nowtri].sides[(nownode + 2)%3], false, false);

			tmptri[1].setnode(_nodenum, this->elements[_nowtri].nodes[(nownode + 1)%3], this->elements[neitri].nodes[neinode]);
			tmptri[1].setneighbor(this->elements[neitri].neighbors[(neinode + 1)%3], this->elements.size() + 1, _nowtri);
			tmptri[1].setside(this->elements[neitri].sides[(neinode + 1)%3], false, false);

			tmptri[2].setnode(_nodenum, this->elements[_nowtri].nodes[(nownode + 2)%3], this->elements[_nowtri].nodes[nownode]);
			tmptri[2].setneighbor(this->elements[_nowtri].neighbors[(nownode + 1)%3], _nowtri, this->elements.size() + 1);
			tmptri[2].sides[0] = this->elements[_nowtri].sides[(nownode + 1)%3];

			tmptri[3].setnode(_nodenum, this->elements[neitri].nodes[neinode], this->elements[_nowtri].nodes[(nownode + 2)%3]);
			tmptri[3].setneighbor(this->elements[neitri].neighbors[(neinode + 2)%3], this->elements.size(), neitri);
			tmptri[3].sides[0] = this->elements[neitri].sides[(neinode + 2)%3];

			int nei1 = this->elements[_nowtri].neighbors[(nownode + 1)%3];
			if (nei1 != -1) {
				this->elements[nei1].neighbors[this->elements[nei1].oppositenode(_nowtri)] = this->elements.size();
			}

			int nei2 = this->elements[neitri].neighbors[(neinode + 2)%3];
			if (nei2 != -1) {
				this->elements[nei2].neighbors[this->elements[nei2].oppositenode(neitri)] = this->elements.size() + 1;
			}

			if (tmptri[0].nodes[1] == _nodenumm1 || tmptri[0].nodes[1] == _nodenump1) {
				tmptri[2].sides[1] = true;
				tmptri[0].sides[2] = true;
			}
			if (tmptri[1].nodes[1] == _nodenumm1 || tmptri[1].nodes[1] == _nodenump1) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
			}
			if (tmptri[2].nodes[1] == _nodenumm1 || tmptri[2].nodes[1] == _nodenump1) {
				tmptri[2].sides[2] = true;
				tmptri[3].sides[1] = true;
			}
			if (tmptri[3].nodes[1] == _nodenumm1 || tmptri[3].nodes[1] == _nodenump1) {
				tmptri[3].sides[2] = true;
				tmptri[1].sides[1] = true;
			}

			if (this->elements[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				tmptri[2].sides[2] = true;
				tmptri[3].sides[1] = true;
				this->nodes[_nodenum].isonboundary = true;
			}

			tmptri[0].getangle(this->nodes);
			tmptri[1].getangle(this->nodes);
			tmptri[2].getangle(this->nodes);
			tmptri[3].getangle(this->nodes);

			stack.push_back(_nowtri);
			stack.push_back(neitri);
			stack.push_back(this->elements.size());
			stack.push_back(this->elements.size() + 1);

			this->elements[_nowtri].copy(tmptri[0]);
			this->elements[neitri].copy(tmptri[1]);
			this->elements.push_back(tmptri[2]);
			this->elements.push_back(tmptri[3]);
		}
		
		//there is not neighbor element
		else {
			std::array<Element<T>, 2> tmptri;

			tmptri[0].setnode(_nodenum, this->elements[_nowtri].nodes[nownode], this->elements[_nowtri].nodes[(nownode + 1)%3]);
			tmptri[0].setneighbor(this->elements[_nowtri].neighbors[(nownode + 2)%3], -1, this->elements.size());
			tmptri[0].setside(this->elements[_nowtri].sides[(nownode + 2)%3], false, false);

			tmptri[1].setnode(_nodenum, this->elements[_nowtri].nodes[(nownode + 2)%3], this->elements[_nowtri].nodes[nownode]);
			tmptri[1].setneighbor(this->elements[_nowtri].neighbors[(nownode + 1)%3], _nowtri, -1);
			tmptri[1].setside(this->elements[_nowtri].sides[(nownode + 1)%3], false, false);

			int nei1 = this->elements[_nowtri].neighbors[(nownode + 1)%3];
			if (nei1 != -1) {
				this->elements[nei1].neighbors[this->elements[nei1].oppositenode(_nowtri)] = this->elements.size();
			}

			if (tmptri[0].nodes[1] == _nodenumm1 || tmptri[0].nodes[1] == _nodenump1) {
				tmptri[1].sides[1] = true;
				tmptri[0].sides[2] = true;
			}
			if (tmptri[0].nodes[2] == _nodenumm1 || tmptri[0].nodes[2] == _nodenump1) {
				tmptri[0].sides[1] = true;
			}
			if (tmptri[1].nodes[1] == _nodenumm1 || tmptri[1].nodes[1] == _nodenump1) {
				tmptri[1].sides[2] = true;
			}

			if (this->elements[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				this->nodes[_nodenum].isonboundary = true;
			}

			tmptri[0].getangle(this->nodes);
			tmptri[1].getangle(this->nodes);

			stack.push_back(_nowtri);
			stack.push_back(this->elements.size());

			this->elements[_nowtri].copy(tmptri[0]);
			this->elements.push_back(tmptri[1]);
		}
		
		this->swapping(stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void Delaunay<T>::swapping(std::vector<int>& _stack, int _nodenump1, int _nodenumm1) {
		while (_stack.size() > 0) {
			int nowstack = _stack[_stack.size() - 1];
			_stack.pop_back();
			
			int neighbortri = this->elements[nowstack].neighbors[0];
			
			if (neighbortri >= 0 && this->elements[neighbortri].active == true) {
				int neighbornode = this->elements[neighbortri].oppositenode(nowstack);
				T r0 = this->nodes[this->elements[nowstack].nodes[1]].distance(this->nodes[this->elements[nowstack].nodes[2]]);
				T r1 = this->nodes[this->elements[nowstack].nodes[0]].distance(this->nodes[this->elements[neighbortri].nodes[neighbornode]]);
				if ((r0 > r1
					&& this->elements[nowstack].inouton(this->elements[neighbortri].nodes[neighbornode], this->nodes) == -1
					&& this->elements[neighbortri].inouton(this->elements[nowstack].nodes[0], this->nodes) == -(neighbornode + 1)
					&& this->elements[nowstack].sides[0] == false)
					|| (this->elements[neighbortri].nodes[neighbornode] == _nodenump1)
					|| (this->elements[neighbortri].nodes[neighbornode] == _nodenumm1)) {
					
					Element<T> tmpelement = Element<T>(this->elements[neighbortri]);

					int neighbor1 = tmpelement.neighbors[(neighbornode + 1)%3];
					if (neighbor1 >= 0) {
						this->elements[neighbor1].neighbors[this->elements[neighbor1].oppositenode(neighbortri)] = nowstack;
					}

					int neighbor2 = this->elements[nowstack].neighbors[1];
					if (neighbor2 >= 0) {
						this->elements[neighbor2].neighbors[this->elements[neighbor2].oppositenode(nowstack)] = neighbortri;
					}

					this->elements[neighbortri].setside(tmpelement.sides[(neighbornode + 2)%3], this->elements[nowstack].sides[1], false);
					this->elements[neighbortri].setnode(this->elements[nowstack].nodes[0], tmpelement.nodes[neighbornode], this->elements[nowstack].nodes[2]);
					this->elements[neighbortri].setneighbor(tmpelement.neighbors[(neighbornode + 2)%3], this->elements[nowstack].neighbors[1], nowstack);

					this->elements[nowstack].setside(tmpelement.sides[(neighbornode + 1)%3], false, this->elements[nowstack].sides[2]);
					this->elements[nowstack].setnode(this->elements[nowstack].nodes[0], this->elements[nowstack].nodes[1], tmpelement.nodes[neighbornode]);
					this->elements[nowstack].setneighbor(tmpelement.neighbors[(neighbornode + 1)%3], neighbortri, this->elements[nowstack].neighbors[2]);

					if (this->elements[nowstack].nodes[2] == _nodenumm1 || this->elements[nowstack].nodes[2] == _nodenump1) {
						this->elements[nowstack].sides[1] = true;
						this->elements[neighbortri].sides[2] = true;
					}

					this->elements[nowstack].getangle(this->nodes);
					this->elements[neighbortri].getangle(this->nodes);

					_stack.push_back(nowstack);
					_stack.push_back(neighbortri);
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::getinternalelement(T _maxside) {
		for (int i = 0; i < this->ADDITIONALNODENUM0; i++) {
			T maxside = this->nodes[this->elements[0].nodes[0]].distance(this->nodes[this->elements[0].nodes[1]]);
			int maxelement = 0, maxnode = 0;
			for (int j = 0; j < this->elements.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (maxside < this->nodes[this->elements[j].nodes[(k + 1)%3]].distance(this->nodes[this->elements[j].nodes[(k + 2)%3]]) && this->elements[j].active == true) {
						maxside = this->nodes[this->elements[j].nodes[(k + 1)%3]].distance(this->nodes[this->elements[j].nodes[(k + 2)%3]]);
						maxelement = j;
						maxnode = k;
					}
				}
			}

			if (maxside < _maxside) {
				break;
			}

			this->nodes.push_back(Node<T>(0.5*(this->nodes[this->elements[maxelement].nodes[(maxnode + 1)%3]].x + this->nodes[this->elements[maxelement].nodes[(maxnode + 2)%3]].x)));

			this->getelementon(maxelement, maxnode, -2, this->nodes.size() - 1, -2);
		}
	}


    template<class T>
    std::vector<Vector<T> > Delaunay<T>::GenerateNodes() {
        std::vector<Vector<T> > x = std::vector<Vector<T> >();
        for (auto node : this->nodes) {
            x.push_back(node.x);
        } 
        return x;
    }
    
    
    template<class T>
    std::vector<std::vector<int> > Delaunay<T>::GenerateElements() {
        std::vector<std::vector<int> > es = std::vector<std::vector<int> >();
        for (auto element : this->elements) {
            es.push_back({ element.nodes[0], element.nodes[1], element.nodes[2] });
        }
        return es;
    }
}