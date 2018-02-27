#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "Matrix.h"
#include "mpi.h"

using Matrix_ns::Matrix;

#define DOUBLE 1
#define FLOAT 2

#define DTYPE DOUBLE

#if DTYPE == DOUBLE
#define dtype double
MPI_Datatype mpi_datatype = MPI_DOUBLE;
#else
#define dtype float
MPI_Datatype mpi_datatype = MPI_FLOAT;
#endif

#define EPS 1e-12


template<class T>
dtype scalar_prod(const T* a, const T* b, size_t size){
    float s = 0;
    for (size_t i = 0; i < size; ++i){
        s += a[i]*b[i];
    }
    return s;
}

template<class T>
dtype norm(const T* a, size_t size) {
    return std::sqrt(scalar_prod(a, a, size));
}


int rank = 0, comm_size;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    try {
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        struct timeval st, et;
        gettimeofday(&st, NULL);

        bool need_rand = false;

        if (argc == 3){
            need_rand = true;
        } else if (argc != 2) {
            throw std::string("Wrong args");
        }

        Matrix<dtype> A;

        std::ifstream in;
        if (!need_rand) {
            in.open(argv[1], std::ios::in);
            A = Matrix<dtype>(in, Matrix_ns::Normal, Matrix_ns::ColMaj);
        } else {
            A = Matrix<dtype>(std::atoi(argv[1]), std::atoi(argv[2]), Matrix_ns::ColMaj, true);
        }

        size_t size = A.n_rows();
        dtype* x_vec = new dtype[size];

        Matrix<dtype> I(size, size, 0.0, Matrix_ns::ColMaj);
        for (size_t i = 0; i < size; ++i) {
            I(i, i) = 1.0;
        }
        Matrix<dtype> U;

        for (size_t ind = 0; ind < A.n_rows(); ++ind) {

            dtype *a_vec = A.get_col(ind);

            std::copy(a_vec, a_vec + size, x_vec);
            dtype a_norm = norm(x_vec+ind, size-ind);

            for (size_t j = 0; j < ind; ++j) {
                x_vec[j] = 0;
            }

            x_vec[ind] -= a_norm;
            dtype x_norm = norm(x_vec, size);

            if (x_norm >= EPS) {
                for (size_t i = 0; i < size; ++i) {
                    x_vec[i] /= x_norm;
                }

                Matrix<dtype> Outer(size, size, 0.0, Matrix_ns::ColMaj);
                for (size_t i = 0; i < size; ++i)
                    for (size_t j = 0; j < size; ++j)
                        Outer(i, j) = 2 * x_vec[i] * x_vec[j];

                U = I - Outer;
            } else {
                U = I;
            }

            A = U*A;
        }

        A.print();

        x_vec[A.n_rows()-1] = - (A(A.n_rows()-1, A.n_cols()-1)/A(A.n_rows()-1, A.n_cols()-2));
        for (int i = (int) A.n_rows()-2; i >= 0; --i){
            dtype sum = A(i, A.n_cols()-1);
            for (int j = i+1; j < A.n_cols(); ++j) {
                sum += A(i, j)*x_vec[j];
            }
            x_vec[i] = - sum/A(i,i);
        }

        gettimeofday(&et, NULL);
        int elapsed = ((et.tv_sec - st.tv_sec) * 1000000) + (et.tv_usec - st.tv_usec);

        for (size_t ind = 0; ind < A.n_rows(); ++ind) {
            std::cout << "x_" << ind << ": " << -x_vec[ind] << std::endl;
        }
        if (rank == 0)
            std::cout << "Time (microsec): " << elapsed << std::endl;
        delete x_vec;
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


