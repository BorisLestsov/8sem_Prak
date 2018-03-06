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

#define EPS 1e-12


dtype f(size_t i, size_t j){
    //return rand() / (dtype) RAND_MAX;
    //return ((1+i+j)%20)/((1+i+j)%15);
    return (i + j) / (i+j+1.0);
}


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

        struct timeval st_1, et_1, st_2, et_2;

        bool need_rand = false;

        if (argc == 3){
            need_rand = true;
        } else if (argc != 2) {
            throw std::string("Wrong args");
        }

        size_t glob_rows,
                glob_cols,
                my_cols = 0;

        std::ifstream in;
        if (!need_rand) {
            in.open(argv[1], std::ios::in);
            in >> glob_rows >> glob_cols;
        } else {
            srand(time(0));
            glob_rows = std::atoi(argv[1]);
            glob_cols = std::atoi(argv[2]);
        }

        my_cols = glob_cols / comm_size;
        if (glob_cols % comm_size != 0){
            if (rank < glob_cols % comm_size)
                my_cols += 1;
        }

        int my_cols_offset;
        Matrix<dtype> A(glob_rows, my_cols, 0.0, Matrix_ns::ColMaj);
        for (size_t i = 0; i < glob_rows; ++i){
            my_cols_offset = rank;
            for (size_t j = 0; j < glob_cols; ++j){
                dtype tmp;
                if (!need_rand) {
                    in >> tmp;
                } else {
                    tmp = f(i, j); 
                }
                if (j % comm_size == rank) {
                    std::cout << rank << "  " << i << "  " << j << "  " << my_cols_offset << "  " << tmp << std::endl;
                    A(i, j - my_cols_offset) = tmp;

                    my_cols_offset += comm_size-1;
                }
            }
        }


        Matrix<dtype> Acopy = A;

        gettimeofday(&st_1, NULL);

        size_t size = glob_rows;
        dtype* x_vec = new dtype[size];

        my_cols_offset = 0;
        for (size_t ind = 0; ind < glob_rows; ++ind) {
            int root = ind % comm_size;
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
                my_cols_offset += comm_size;
            }
            MPI_Bcast(x_vec, size, mpi_datatype, root, MPI_COMM_WORLD);
            for (size_t col_i = 0; col_i < my_cols; ++col_i) {
                dtype *y_vec = A.get_col(col_i);
                dtype scal = 2*scalar_prod(x_vec, y_vec, size);
                for (size_t j = 0; j < size; ++j) {
                    y_vec[j] -= scal*x_vec[j];
                }
            }
        }

        gettimeofday(&et_1, NULL);

        gettimeofday(&st_2, NULL);

        dtype* tmp_vec = new dtype[size+1];
        dtype* b_vec = new dtype[size];

        if (rank == (glob_cols - 1) % comm_size){
            dtype* ptr = A.get_col(A.n_cols() - 1);
            std::copy(ptr, ptr + size, b_vec);
        }
        MPI_Bcast(b_vec, size, mpi_datatype, (glob_cols - 1) % comm_size, MPI_COMM_WORLD);

        my_cols_offset = 0;
        for (int ind = size-1; ind >= 0; --ind) {
            dtype val;
            dtype* ptr = tmp_vec;
            int root = ind % comm_size;
            if (rank == root) {
                val = b_vec[ind] / A(ind, ind - my_cols_offset);
                tmp_vec = A.get_col(ind - my_cols_offset);
            }
            MPI_Bcast(&val, 1, mpi_datatype, root, MPI_COMM_WORLD);
            MPI_Bcast(tmp_vec, ind, mpi_datatype, root, MPI_COMM_WORLD);

            x_vec[ind] = val;

            for (int loc_row = rank; loc_row < ind; loc_row += comm_size) {
                b_vec[loc_row] -= val*tmp_vec[loc_row];
            }

            tmp_vec = ptr;

            my_cols_offset += comm_size;
        }
        

        gettimeofday(&et_2, NULL);
/*
        dtype* b_vec, *res_vec = new dtype[glob_rows];
        if (rank == glob_cols % comm_size)
            b_vec = Acopy.get_col(Acopy.n_cols()-1);

 
        if (rank == comm_size-1)
            my_cols -= 1;

        for (int i = 0; i < size; ++i){
            res_vec[i] = 0.0;
            tmp_vec[i] = 0.0;
            for (int j = 0; j < my_cols; ++j){
                tmp_vec[i] += Acopy(i, j) * x_vec[my_cols_offset+j];
            }
        }

        MPI_Reduce(tmp_vec, res_vec, size, mpi_datatype, MPI_SUM, comm_size-1, MPI_COMM_WORLD);
       
        if (rank == comm_size-1)
            for (size_t i = 0; i < size; ++i){
                res_vec[i] -= b_vec[i];
            }
*/
        dtype diff = 10;
  //      if (rank == comm_size - 1)
  //          diff = norm(res_vec, size);
        
        int elapsed_1 = ((et_1.tv_sec - st_1.tv_sec) * 1000000) + (et_1.tv_usec - st_1.tv_usec);
        int elapsed_2 = ((et_2.tv_sec - st_2.tv_sec) * 1000000) + (et_2.tv_usec - st_2.tv_usec);

        if (rank == comm_size-1)
            for (size_t ind = 0; ind < size; ++ind) {
                std::cout << "x_" << ind << ": " << x_vec[ind] << std::endl;
           }
        if (rank == comm_size-1)
            std::cout <<  "Mat_size " << size << " Comm_size " << comm_size << " Forward_Time_(microsec) " << elapsed_1 << "  Backward_Time_(microsec) " << elapsed_2 << " diff " << diff << std::endl;
        delete[] x_vec;
        delete[] tmp_vec;

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
