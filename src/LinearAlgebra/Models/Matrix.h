//*****************************************************************************
//  Title       :src/LinearAlgebra/Matrix.h
//  Author      :Tanabe Yuta
//  Date        :2019/11/07
//  Copyright   :(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <iostream>
#include <cassert>


#include "Vector.h"


namespace PANSFEM2{
    template<class T>
    class Vector;


    template<class T>
    class Matrix{
public:
        Matrix();
        virtual ~Matrix();
        Matrix(int _row, int _col);
        Matrix(const Matrix<T>& _mat);


        int ROW() const;
        int COL() const;
        T& operator()(int _i, int _j);


        Matrix<T>& operator=(const Matrix<T>& _mat);
        Matrix<T>& operator+=(const Matrix<T>& _mat);
        Matrix<T>& operator-=(const Matrix<T>& _mat);
        Matrix<T>& operator*=(T _a);
        Matrix<T>& operator/=(T _a);


        Matrix<T> operator+(const Matrix<T>& _mat);
        Matrix<T> operator-(const Matrix<T>& _mat);
        Matrix<T> operator-();
        Matrix<T> operator*(const Matrix<T>& _mat);
        Matrix<T> operator*(T _a);
        Matrix<T> operator/(T _a);


        template<class U>
	    friend std::ostream& operator << (std::ostream& _out, const Matrix<U>& _mat);
        template<class U>
        friend Vector<U> operator*(const Matrix<U>& _mat, const Vector<U>& _vec);


        Matrix<T> Transpose();
        T Determinant();
        Matrix<T> Inverse();
        Matrix<T> Cofactor(int _i, int _j);
        Matrix<T> Vstack(const Matrix<T>& _mat);
        Matrix<T> Hstack(const Matrix<T>& _mat);


        template<class U>
        friend Matrix<U> Vector<U>::Transpose();
        template<class U>
        friend Matrix<U> Vector<U>::operator*(const Matrix<U>& _mat);


        template<class U>
        friend Matrix<U> Identity(int _row);


protected:
        int row;            //Number of row
        int col;            //Number of column
        T* values;          //Values of matrix
    };


    template<class T>
    Matrix<T>::Matrix() {
        this->row = 0;
        this->col = 0;
    }


    template<class T>
    Matrix<T>::~Matrix(){
        delete[] this->values;
    }


    template<class T>
    Matrix<T>::Matrix(int _row, int _col) {
        this->row = _row;
        this->col = _col;
        this->values = new T[this->row * this->col]();
    }


    template<class T>
    Matrix<T>::Matrix(const Matrix<T>& _mat) {
        this->row = _mat.row;
        this->col = _mat.col;
        this->values = new T[this->row * this->col];
        for(int i = 0; i < this->row * this->col; i++){
            this->values[i] = _mat.values[i];
        }
    }


    template<class T>
    int Matrix<T>::ROW() const {
        return this->row;
    }


    template<class T>
    int Matrix<T>::COL() const {
        return this->col;
    }


    template<class T>
    T& Matrix<T>::operator()(int _i, int _j){
        assert(0 <= _i && _i < this->row && 0 <= _j && _j < this->col);
        return this->values[_i * this->col + _j];
    }


    template<class T>
    Matrix<T>& Matrix<T>::operator=(const Matrix<T>& _mat){
        if(this != &_mat){
            if(this->row != 0 && this->col != 0){
                delete[] this->values;
            }

            this->row = _mat.row;
            this->col = _mat.col;
            this->values = new T[this->row * this->col];
            for(int i = 0; i < this->row * this->col; i++){
                this->values[i] = _mat.values[i];
            }
        }
        return *this;
    }


    template<class T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& _mat){
        assert(this->row == _mat.row && this->col == _mat.col);
        for(int i = 0; i < this->row * this->col; i++){
            this->values[i] += _mat.values[i];
        }
        return *this;
    }


    template<class T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& _mat){
        assert(this->row == _mat.row && this->col == _mat.col);
        for(int i = 0; i < this->row * this->col; i++){
            this->values[i] -= _mat.values[i];
        }
        return *this;
    }


    template<class T>
    Matrix<T>& Matrix<T>::operator*=(T _a){
        for(int i = 0; i < this->row * this->col; i++){
            this->values[i] *= _a;
        }
        return *this;
    }


    template<class T>
    Matrix<T>& Matrix<T>::operator/=(T _a){
        for(int i = 0; i < this->row * this->col; i++){
            this->values[i] /= _a;
        }
        return *this;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator+(const Matrix<T>& _mat){
        assert(this->row == _mat.row && this->col == _mat.col);
        Matrix<T> mat = *this;
        for(int i = 0; i < mat.row * mat.col; i++){
            mat.values[i] += _mat.values[i];
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator-(const Matrix<T>& _mat){
        assert(this->row == _mat.row && this->col == _mat.col);
        Matrix<T> mat = *this;
        for(int i = 0; i < mat.row * mat.col; i++){
            mat.values[i] -= _mat.values[i];
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator-(){
        Matrix<T> mat = *this;
        for(int i = 0; i < mat.row * mat.col; i++){
            mat.values[i] *= -1.0;
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator*(const Matrix<T>& _mat){
        assert(this->col == _mat.row);
        Matrix<T> mat = Matrix<T>(this->row, _mat.col);
        for(int i = 0; i < mat.row; i++){
            for(int j = 0; j < mat.col; j++){
                mat.values[i * mat.col + j] = T();
                for(int k = 0; k < this->col; k++){
                    mat.values[i * mat.col + j] += this->values[i * this->col + k] * _mat.values[k * _mat.col + j];
                }
            }
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator*(T _a){
        Matrix<T> mat = *this;
        for(int i = 0; i < mat.row * mat.col; i++){
            mat.values[i] *= _a;
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::operator/(T _a){
        Matrix<T> mat = *this;
        for(int i = 0; i < mat.row * mat.col; i++){
            mat.values[i] /= _a;
        }
        return mat;
    }


    template<class U>
    std::ostream& operator << (std::ostream& _out, const Matrix<U>& _mat){
        for(int i = 0; i < _mat.row; i++){
            for(int j = 0; j < _mat.col; j++){
                _out << _mat.values[i * _mat.col + j] << "\t";
            }
            _out << std::endl;
        }
        return _out;
    }


    template<class U>
    Vector<U> operator*(const Matrix<U>& _mat, const Vector<U>& _vec){
        assert(_mat.col == _vec.size);
        Vector<U> vec = Vector<U>(_mat.row);
        for(int i = 0; i < vec.size; i++){
            for(int j = 0; j < _mat.col; j++){
                vec(i) += _mat.values[i * _mat.col + j] * _vec.values[j];
            }
        }
        return vec;
    }


    template<class T>
    Matrix<T> Matrix<T>::Transpose(){
        Matrix<T> mat = Matrix<T>(this->col, this->row);
        for(int i = 0; i < mat.row; i++){
            for(int j = 0; j < mat.col; j++){
                mat.values[i * mat.col + j] =  this->values[j * this->col + i];
            }
        }
        return mat;
    }


    template<class T>
    T Matrix<T>::Determinant(){
        assert(this->row == this->col && this->row != 0);
        if(this->row == 1) {
            return this->values[0];
        } else if(this->row == 2) {
            return this->values[0]*this->values[3] - this->values[1]*this->values[2];
        } else if(this->row == 3) {
            return - this->values[8]*this->values[1]*this->values[3] - this->values[7]*this->values[5]*this->values[0] - this->values[2]*this->values[4]*this->values[6]
                    + this->values[6]*this->values[1]*this->values[5] + this->values[7]*this->values[3]*this->values[2] + this->values[0]*this->values[4]*this->values[8];
        } else {
            T value = T();
            for(int i = 0; i < this->row; i++){
                value += pow(-1.0, i) * this->values[i * this->col] * this->Cofactor(i, 0).Determinant();
            }
            return value;
        }
    }


    template<class T>
    Matrix<T> Matrix<T>::Inverse(){
        assert(this->row == this->col);
        Matrix<T> mat = Matrix<T>(this->row, this->col);
        for(int i = 0; i < this->row; i++){
            for(int j = 0; j < this->col; j++){
                mat.values[i * mat.col + j] = pow(-1.0, i + j) * this->Cofactor(j, i).Determinant();
            }
        }
        return mat / this->Determinant();
    }


    template<class T>
    Matrix<T> Matrix<T>::Cofactor(int _i, int _j){
        assert(0 <= _i && _i < this->row && 0 <= _j && _j < this->col);
        Matrix<T> mat = Matrix<T>(this->row - 1, this->col - 1);
        for(int i = 0; i < mat.row; i++){
            for(int j = 0; j < mat.col; j++){
                if(i < _i){
                    if(j < _j){
                        mat.values[i * mat.col + j] = this->values[i * this->col + j];
                    } else {
                        mat.values[i * mat.col + j] = this->values[i * this->col + (j + 1)];
                    }
                } else {
                    if(j < _j){
                        mat.values[i * mat.col + j] = this->values[(i + 1) * this->col + j];
                    } else {
                        mat.values[i * mat.col + j] = this->values[(i + 1) * this->col + (j + 1)];
                    }
                }
            }
        }
        return mat;
    }


    template<class T>
    Matrix<T> Matrix<T>::Vstack(const Matrix<T>& _mat){
        assert(this->col == _mat.col);
        Matrix<T> mat = Matrix<T>(this->row + _mat.row, this->col);
        for(int i = 0; i < this->row; i++){
            for(int j = 0; j < this->col; j++){
                mat.values[i * mat.col + j] = this->values[i * this->col + j];
            }
        }
        for(int i = 0; i < _mat.row; i++){
            for(int j = 0; j < _mat.col; j++){
                mat.values[(i + this->row) * mat.col + j] = _mat.values[i * _mat.col + j];
            }
        }
        return mat;
    }
    
    
    template<class T>
    Matrix<T> Matrix<T>::Hstack(const Matrix<T>& _mat){
        assert(this->row == _mat.row);
        Matrix<T> mat = Matrix<T>(this->row, this->col + _mat.col);
        for(int i = 0; i < this->row; i++){
            for(int j = 0; j < this->col; j++){
                mat.values[i * mat.col + j] = this->values[i * this->col + j];
            }
        }
        for(int i = 0; i < _mat.row; i++){
            for(int j = 0; j < _mat.col; j++){
                mat.values[i * mat.col + (j + this->col)] = _mat.values[i * _mat.col + j];
            }
        }
        return mat;
    }


    template<class U>
    Matrix<U> Identity(int _row){
        Matrix<U> mat = Matrix<U>(_row, _row);
        for(int i = 0; i < mat.row; i++){
            mat.values[i * mat.col + i] = 1.0;
        }
        return mat;
    }
}