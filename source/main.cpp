#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "Matrix.h"
#include "mpi.h"

using Matrix_ns::Matrix;

#define DOUBLE 1
#define FLOAT 2

#define DTYPE FLOAT

#if DTYPE == DOUBLE
#define dtype double
MPI_Datatype mpi_datatype = MPI_DOUBLE;
#else
#define dtype float
MPI_Datatype mpi_datatype = MPI_FLOAT;
#endif

#define EPS 1e-8


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

        size_t glob_rows,
                glob_cols,
                my_cols,
                cols_stride,
                my_cols_offset,
                next_offset;

        std::ifstream in;
        if (!need_rand) {
            in.open(argv[1], std::ios::in);
            in >> glob_rows >> glob_cols;
        } else {
            srand(time(0));
            glob_rows = std::atoi(argv[1]);
            glob_cols = std::atoi(argv[2]);
        }

        cols_stride = glob_cols/comm_size;
        my_cols = cols_stride;
        if (rank == comm_size-1)
            my_cols += glob_cols % comm_size;

        my_cols_offset = rank * cols_stride;
        next_offset = my_cols_offset+my_cols;

//        if (rank == 3){
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
                if (!need_rand) {
                    in >> tmp;
                } else {
                    tmp = rand() / (dtype) RAND_MAX;
                }
                if (j >= my_cols_offset && j < next_offset) {
                    A(i, j - my_cols_offset) = tmp;
                }
            }
        }

        
        size_t size = glob_rows;
        dtype* x_vec = new dtype[size];

        for (size_t ind = 0; ind < glob_rows; ++ind) {
            int root = std::min(ind / cols_stride, (size_t) comm_size-1);
            if (rank == root) {

                dtype *a_vec = A.get_col(ind - my_cols_offset);

                std::copy(a_vec, a_vec + size, x_vec);
                dtype a_norm = norm(x_vec + ind, size - ind);


                for (size_t j = 0; j < ind; ++j) {
                    x_vec[j] = 0;
                }

                x_vec[ind] -= a_norm;
                dtype x_norm = norm(x_vec, size);

                if (x_norm >= EPS) {
                    for (size_t i = 0; i < size; ++i) {
                        x_vec[i] /= x_norm;
                    }
                }
            }
            MPI_Bcast(x_vec, size, mpi_datatype, root, MPI_COMM_WORLD);
            for (int col_i = 0; col_i < my_cols; ++col_i) {
                dtype *y_vec = A.get_col(col_i);
                dtype scal = 2*scalar_prod(x_vec, y_vec, size);
                for (int j = 0; j < size; ++j) {
                    y_vec[j] -= scal*x_vec[j];
                }
            }
        }

        Matrix<dtype> Afin;
        if (rank == 0)
            Afin = Matrix<dtype>(glob_rows, glob_cols, 0.0, Matrix_ns::ColMaj);

        MPI_Gather(A.data(), glob_rows*cols_stride, mpi_datatype, Afin.data(), glob_rows*cols_stride, mpi_datatype, 0, MPI_COMM_WORLD);

        if (comm_size != 1) {
            if (rank == comm_size - 1)
                MPI_Send(A.data() + cols_stride * glob_rows,
                         glob_rows * (glob_cols - cols_stride * comm_size),
                         mpi_datatype, 0, 0, MPI_COMM_WORLD);
            if (rank == 0)
                MPI_Recv(Afin.data() + (cols_stride * comm_size) * glob_rows,
                         glob_rows * (glob_cols - cols_stride * comm_size),
                         mpi_datatype, comm_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
//        if (rank==0)
//            Afin.print();


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
        }
        
        gettimeofday(&et, NULL);
        int elapsed = ((et.tv_sec - st.tv_sec) * 1000000) + (et.tv_usec - st.tv_usec);

        //for (size_t ind = 0; ind < Afin.n_rows(); ++ind) {
        //    std::cout << "x_" << ind << ": " << -x_vec[ind] << std::endl;
        //}
        if (rank == 0)
            std::cout << "Size: " << comm_size << "; Time (microsec): " << elapsed << std::endl;
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


