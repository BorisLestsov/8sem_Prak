#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <complex>

typedef unsigned int uint;

namespace Matrix_ns {

enum OutType {O_FXD, O_SCI};
enum InType {Normal, Binary};
enum ElOrder {RowMaj=0, ColMaj=1};

using namespace std;


template<class T>
class Matrix {
    T *arr;
    size_t rows;
    size_t cols;

    ElOrder ord;

public:
    Matrix();
    Matrix(T* buf, size_t rows, size_t cols, ElOrder order = RowMaj);
    Matrix(size_t rows, size_t cols, bool is_rand = false, ElOrder order = RowMaj);
    Matrix(size_t rows, size_t cols, T elem, ElOrder order = RowMaj);
    Matrix(const Matrix<T> &matr);
    //Matrix(Matrix&& m);

    Matrix(istream& in, InType in_t = Normal, ElOrder order = RowMaj);

    ~Matrix();


    inline size_t n_rows() const ;
    inline size_t n_cols() const ;
    inline ElOrder order() const ;
    inline size_t size() const;



    Matrix<T> operator+(const Matrix &matr) throw(string);
    Matrix<T> operator-(const Matrix &matr) throw(string);
    Matrix<T> operator*(const Matrix &matr) throw(string);

    inline T &operator()(size_t i, size_t j) const throw(string);

    Matrix<T>& operator=(const Matrix &m);
    //Matrix<T>& operator=(Matrix&& m);

    T* data() const;

    void write(ostream& o) const;
    void print(OutType o_type = O_FXD,
               uint precision = 8,
               ostream &stream = cout) const throw(string);

    double norm();

    template<class Y>
    friend const istream &operator>>(const istream &strm, Matrix<Y> &matr);
};


template<class T>
class Matrix< std::complex<T> > {
    std::complex<T> *arr;
    size_t rows;
    size_t cols;

public:
    Matrix();
    Matrix(std::complex<T>* buf, size_t rows, size_t cols);
    Matrix(size_t rows, size_t cols, bool is_rand = false);
    Matrix(size_t rows, size_t cols, std::complex<T> elem);
    Matrix(const Matrix<std::complex<T> > &matr);
    //Matrix(Matrix&& m);

    Matrix(istream& in, InType in_t = Normal);

    ~Matrix();


    inline size_t n_rows() const ;
    inline size_t n_cols() const ;
    inline size_t size() const;

    Matrix<std::complex<T> > operator+(const Matrix &matr) throw(string);
    Matrix<std::complex<T> > operator-(const Matrix &matr) throw(string);
    Matrix<std::complex<T> > operator*(const Matrix &matr) throw(string);

    inline std::complex<T> &operator()(size_t i, size_t j) const throw(string);

    Matrix<std::complex<T> >& operator=(const Matrix &m);
    //Matrix<T>& operator=(Matrix&& m);

    std::complex<T>* data() const;

    void write(ostream& o) const;
    void print(OutType o_type = O_FXD,
               uint precision = 8,
               ostream &stream = cout) const throw(string);

    double norm();

    template<class Y>
    friend const istream &operator>>(const istream &strm, Matrix<std::complex<Y> > &matr);
};


    //-----implementation-------
#include "Matrix.hpp"

}

#endif // MATRIX_H_INCLUDED
