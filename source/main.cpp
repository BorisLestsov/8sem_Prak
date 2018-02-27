#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

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

        if (argc != 2){
            throw std::string("Wrong args");
        }

        std::ifstream in(argv[1], std::ios::in);

        size_t glob_rows,
                glob_cols,
                my_cols,
                cols_stride,
                my_cols_offset,
                next_offset;

        in >> glob_rows >> glob_cols;
        cols_stride = glob_cols/comm_size;

        my_cols = cols_stride;
        if (rank == comm_size-1)
            my_cols += glob_cols % comm_size;

        my_cols_offset = rank * cols_stride;
        next_offset = my_cols_offset+my_cols;

//        if (rank == 1){
//            std::cout << glob_rows << "  "
//                    << glob_cols << "  "
//                    << my_cols << "  "
//                    << cols_stride << "  "
//                    << my_cols_offset << "  "
//                    << next_offset << std::endl;
//        }


        Matrix<dtype> A(glob_rows, my_cols, 0.0, Matrix_ns::ColMaj);
        for (size_t i = 0; i < glob_rows; ++i){
            for (size_t j = 0; j < glob_cols; ++j){
                dtype tmp;
                in >> tmp;
                if (j >= my_cols_offset && j < next_offset) {
                    A(i, j - my_cols_offset) = tmp;
                }
            }
        }

        size_t size = glob_rows;
        dtype* x_vec = new dtype[size];

        Matrix<dtype> U(glob_rows, glob_rows, 0.0, Matrix_ns::ColMaj);
        for (size_t ind = 0; ind < glob_rows; ++ind) {
            int root = ind / cols_stride;
            //ind >= my_cols_offset && ind < next_offset
            if (rank == root) {
                Matrix<dtype> I(size, size, 0.0, Matrix_ns::ColMaj);
                Matrix<dtype> Outer(size, size, 0.0, Matrix_ns::ColMaj);

                dtype *a_vec = A.get_col(ind - my_cols_offset);

                std::copy(a_vec, a_vec + size, x_vec);
                dtype a_norm = norm(x_vec + ind, size - ind);

                for (size_t j = 0; j < ind; ++j) {
                    x_vec[j] = 0;
                }

                x_vec[ind] -= a_norm;
                dtype x_norm = norm(x_vec, size);

                for (size_t i = 0; i < size; ++i) {
                    x_vec[i] /= x_norm;
                }

                for (size_t i = 0; i < size; ++i) {
                    I(i, i) = 1.0;
                }

                for (size_t i = 0; i < size; ++i)
                    for (size_t j = 0; j < size; ++j)
                        Outer(i, j) = 2 * x_vec[i] * x_vec[j];

                U = I - Outer;
            }

            MPI_Bcast(U.data(), U.n_rows()*U.n_cols(), mpi_datatype, root, MPI_COMM_WORLD);
            A = U * A;
        }

        Matrix<dtype> Afin;
        if (rank == 0)
            Afin = Matrix<dtype>(glob_rows, glob_cols, 0.0, Matrix_ns::ColMaj);

        MPI_Gather(A.data(), glob_rows*cols_stride, mpi_datatype, Afin.data(), glob_rows*cols_stride, mpi_datatype, 0, MPI_COMM_WORLD);

        if (rank==comm_size-1)
            MPI_Send(A.data()+cols_stride*glob_rows,
                     glob_rows*(glob_cols-cols_stride*comm_size),
                     mpi_datatype, 0, 0, MPI_COMM_WORLD);
        if (rank==0)
            MPI_Recv(Afin.data()+(cols_stride*comm_size)*glob_rows,
                     glob_rows*(glob_cols-cols_stride*comm_size),
                     mpi_datatype, comm_size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        if (rank == 0) {
            x_vec[Afin.n_rows() - 1] = -(Afin(Afin.n_rows() - 1, Afin.n_cols() - 1) /
                                        Afin(Afin.n_rows() - 1, Afin.n_cols() - 2));
            for (int i = (int) Afin.n_rows() - 2; i >= 0; --i) {
                dtype sum = Afin(i, Afin.n_cols() - 1);
                for (int j = i + 1; j < Afin.n_cols(); ++j) {
                    sum += Afin(i, j) * x_vec[j];
                }
                x_vec[i] = -sum / Afin(i, i);
            }
            for (size_t ind = 0; ind < A.n_rows(); ++ind) {
                std::cout << "x_" << ind << ": " << -x_vec[ind] << std::endl;
            }
        }

        delete x_vec;

    }
    catch (const std::string& e) {
        std::cerr << "rank " << rank << ": " << e << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Finalize();
    return 0;
}


