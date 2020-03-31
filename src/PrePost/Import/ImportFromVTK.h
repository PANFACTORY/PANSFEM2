//*****************************************************************************
//  Title		:   src/PrePost/Import/ImportFromVTK.h
//  Author  	:   Tanabe Yuta
//  Date		:   2020/03/28
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma oncce
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>


#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //********************Import Model from VTK********************
    template<class T>
    class ImportModelFromVTK{
public:
        ImportModelFromVTK(std::string _fname, int _dimension);
        ~ImportModelFromVTK();
        std::vector<Vector<T> > ImportPOINTS();
        std::vector<std::vector<int> > ImportCELLS();
        std::vector<Vector<T> > ImportPOINTVECTORS(std::string _keyword);
        std::vector<T> ImportPOINTSCALARS(std::string _keyword);
        std::vector<Vector<T> > ImportCELLVECTORS(std::string _keyword);
        std::vector<T> ImportCELLSCALARS(std::string _keyword);


private:
        std::ifstream ifs;      //  File stream
        int dimension;
        int mode = 0;           //  1:import POINT_DATA     2:import CELL_DATA
        int pointssize = 0;
        int cellssize = 0;
    };


    template<class T>
    ImportModelFromVTK<T>::ImportModelFromVTK(std::string _fname, int _dimension){
        assert(0 < _dimension && _dimension < 4);
        this->ifs.open(_fname);
        if (!this->ifs.is_open()) {
			std::cout << "VTK file " << _fname << " open error!" << std::endl;
		}
        this->dimension = _dimension;
    }


    template<class T>
    ImportModelFromVTK<T>::~ImportModelFromVTK(){
        if (!this->ifs.is_open()) {
			this->ifs.close();
		}
    }


    template<class T>
    std::vector<Vector<T> > ImportModelFromVTK<T>::ImportPOINTS(){
        std::string linestrings;
        std::vector<Vector<T> > nodes;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(linestrings.find("POINTS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                nodes = std::vector<Vector<T> >(this->pointssize, Vector<T>(this->dimension));
                for(auto& node : nodes){
                    getline(this->ifs, linestrings);
                    std::stringstream valuestream(linestrings);
                    for(int i = 0; i < this->dimension; i++){
                        valuestream >> node(i);
                    }
                }
                return nodes;
            }
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return nodes;
    }


    template<class T>
    std::vector<std::vector<int> > ImportModelFromVTK<T>::ImportCELLS(){
        std::string linestrings;
        std::vector<std::vector<int> > elements;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(linestrings.find("CELLS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                elements = std::vector<std::vector<int> >(cellssize);
                for(auto& element : elements){
                    getline(this->ifs, linestrings);
                    std::stringstream valuestream(linestrings);
                    int elementsize;
                    valuestream >> elementsize;
                    element = std::vector<int>(elementsize);
                    for(auto& node : element){
                        valuestream >> node;
                    }
                }
                return elements;
            }    
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return elements;
    }


    template<class T>
    std::vector<Vector<T> > ImportModelFromVTK<T>::ImportPOINTVECTORS(std::string _keyword){
        std::string linestrings;
        std::vector<Vector<T> > pointvectors;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(this->mode == 1 && linestrings.find("VECTORS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> linestrings;
                if(linestrings == _keyword){
                    pointvectors = std::vector<Vector<T> >(this->pointssize, Vector<T>(this->dimension));
                    for(auto& pointvector : pointvectors){
                        getline(this->ifs, linestrings);
                        std::stringstream valuestream(linestrings);
                        for(int i = 0; i < this->dimension; i++){
                            valuestream >> pointvector(i);
                        }
                    }
                    return pointvectors;
                }
            }
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return pointvectors;
    }


    template<class T>
    std::vector<T> ImportModelFromVTK<T>::ImportPOINTSCALARS(std::string _keyword){
        std::string linestrings;
        std::vector<T> pointscalars;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(this->mode == 1 && linestrings.find("SCALARS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> linestrings;
                if(linestrings == _keyword){
                    pointscalars = std::vector<T>(this->pointssize);
                    getline(this->ifs, linestrings);
                    for(auto& pointscalar : pointscalars){
                        getline(this->ifs, linestrings);
                        std::stringstream valuestream(linestrings);
                        valuestream >> pointscalar;
                    }
                    return pointscalars;
                }
            }
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return pointscalars;
    }


    template<class T>
    std::vector<Vector<T> > ImportModelFromVTK<T>::ImportCELLVECTORS(std::string _keyword){
        std::string linestrings;
        std::vector<Vector<T> > cellvectors;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(this->mode == 2 && linestrings.find("VECTORS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> linestrings;
                if(linestrings == _keyword){
                    cellvectors = std::vector<Vector<T> >(this->cellssize, Vector<T>(this->dimension));
                    for(auto& cellvector : cellvectors){
                        getline(this->ifs, linestrings);
                        std::stringstream valuestream(linestrings);
                        for(int i = 0; i < this->dimension; i++){
                            valuestream >> cellvector(i);
                        }
                    }
                    return cellvectors;
                }
            }
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return cellvectors;
    }


    template<class T>
    std::vector<T> ImportModelFromVTK<T>::ImportCELLSCALARS(std::string _keyword){
        std::string linestrings;
        std::vector<T> cellscalars;
        while (getline(this->ifs, linestrings)) {
            if(linestrings.find("POINT_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->pointssize;
                this->mode = 1;
            } else if(linestrings.find("CELL_DATA") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> this->cellssize;
                this->mode = 2;
            } else if(this->mode == 2 && linestrings.find("SCALARS") != std::string::npos){
                std::stringstream keywordstream(linestrings);
                keywordstream >> linestrings >> linestrings;
                if(linestrings == _keyword){
                    cellscalars = std::vector<T>(this->cellssize);
                    getline(this->ifs, linestrings);
                    for(auto& cellscalar : cellscalars){
                        getline(this->ifs, linestrings);
                        std::stringstream valuestream(linestrings);
                        valuestream >> cellscalar;
                    }
                    return cellscalars;
                }
            }
		}
        this->ifs.clear();
        this->ifs.seekg(0, std::ios_base::beg); 
        return cellscalars;
    }
}