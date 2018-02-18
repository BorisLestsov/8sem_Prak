#include <iostream>
#include <fstream>
#include "Matrix.h"
#include <algorithm>
#include <vector>
#include <string>
#include "mpi.h"


using Matrix_ns::Matrix;

#define DEBUG

template<class T>
float scalar_prod(const T* a, const T* b, size_t size){
    float s = 0;
    for (size_t i = 0; i < size; ++i){
        s += a[i]*b[i];
    }
    return s;
}

template<class T>
float norm(const T* a, size_t size) {
    return std::sqrt(scalar_prod(a, a, size));
}


int rank = 0, comm_size;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    try {
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (argc != 2){
            throw std::string("Wrong parameters");
        }

        std::ifstream in(argv[1], std::ios::in);

        Matrix<float> A(in, Matrix_ns::Normal, Matrix_ns::ColMaj);
        A.print();

        size_t size = A.n_rows();
        float* e_vec = new float[size];
        float* x_vec = new float[size];

        Matrix<float> I(size, size, 0.0);
        Matrix<float> Outer(size, size, 0.0);


        size_t ind = 0;


        float* a_vec = A.get_col(ind);
        e_vec[ind] = 1;


        float a_norm = norm(a_vec, size);
        std::copy(a_vec, a_vec+size, x_vec);
        x_vec[ind] -= a_norm;
        float x_norm = norm(x_vec, size);

        for (size_t i = 0; i < size; ++i){
            x_vec[i] /= x_norm;
        }

        for (size_t i = 0; i < size; ++i){
            I(i, i) = 1.0;
        }

        for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < size; ++j)
                Outer(i, j) = 2 * x_vec[i] * x_vec[j];

        A = (I-Outer)*A;

        A.print();


        delete e_vec;

        //std::vector<float> a0 = A.get_col(0);

        //for (size_t i = 0; i < A.n_rows(); ++i){
        //}


    }
    catch (const std::string& e) {
        if (rank == 0)
            std::cerr << e << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Finalize();
    return 0;
}


