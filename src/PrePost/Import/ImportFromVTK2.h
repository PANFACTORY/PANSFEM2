//*****************************************************************************
//  Title		:   src/PrePost/Import/ImportFromVTK2.h
//  Author  	:   Tanabe Yuta
//  Date		:   2021/01/20
//  Copyright	:   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <iostream>
#include <string>
#include <fstream>

#include "../../LinearAlgebra/Models/Vector.h"

namespace PANSFEM2 {
    template<class T>
    class ImportModelFromVTK {
public:
        ImportModelFromVTK() = delete;
        ImportModelFromVTK(std::string _fname) {
            this->ifs.open(_fname);
        }
        ~ImportModelFromVTK() {}

        std::vector<Vector<T> > GenerateNodes();
        std::vector<std::vector<int> > GenerateElements();

private:
        std::ifstream ifs;
    };

    template<class T>
    std::vector<Vector<T> > ImportModelFromVTK<T>::GenerateNodes() {
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        std::vector<Vector<T> > nodes;
        std::string tmp;
        while (ifs >> tmp) {
            if (tmp == "POINTS") {
                int nodesnum;
                ifs >> nodesnum;
                ifs >> tmp;
                for (int i = 0; i < nodesnum; i++) {
                    T x, y, z;
                    ifs >> x >> y >> z;
                    nodes.push_back({ x, y, z });
                }
            }
        }
        return nodes;
    }

    template<class T>
    std::vector<std::vector<int> > ImportModelFromVTK<T>::GenerateElements() {
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        std::vector<std::vector<int> > elements;
        std::string tmp;
        while (ifs >> tmp) {
            if (tmp == "CELLS") {
                int elementsnum;
                ifs >> elementsnum;
                ifs >> tmp;
                for (int i = 0; i < elementsnum; i++) {
                    int elementnum;
                    ifs >> elementnum;
                    std::vector<int> element;
                    for (int j = 0; j < elementnum; j++) {
                        int nodeid;
                        ifs >> nodeid;
                        element.push_back(nodeid);
                    }
                    elements.push_back(element);
                }
            }
        }
        return elements;
    }
}