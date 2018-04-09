#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <cfloat>

#include "Matrix.h"
#include "mpi.h"

using Matrix_ns::Matrix;

#define DOUBLE 1
#define FLOAT 2

#define DTYPE FLOAT

#if DTYPE == FLOAT
#define dtype double
MPI_Datatype mpi_datatype = MPI_DOUBLE;
#define EPS 10*DBL_MIN
#else
#define dtype float
MPI_Datatype mpi_datatype = MPI_FLOAT;
#define EPS 10*FLT_MIN
#endif


dtype f(size_t i, size_t j){
    //return std::max(i, j);
    return rand() / (dtype) RAND_MAX;
    //return (i + j) / (i+j+1.0);
   // return 1.0;///(i+j+1);
    //return i+j;
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

        double st_1, et_1, st_2, et_2;

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
            if (glob_cols != glob_rows+1)
                throw std::string("wrong args");
        }

        my_cols = glob_cols / comm_size;
        if (glob_cols % comm_size != 0){
            if (rank < glob_cols % comm_size)
                my_cols += 1;
        }

//        dtype* tmp1 = new dtype[glob_rows];
//        dtype* sum_bvec = new dtype[glob_rows];

        int my_cols_offset;
        Matrix<dtype> A(glob_rows, my_cols, 0.0, Matrix_ns::ColMaj);
        for (size_t i = 0; i < glob_rows; ++i){
//
//            tmp1[i] = 0;
//            sum_bvec[i] = 0;


            my_cols_offset = rank;
            for (size_t j = 0; j < glob_cols; ++j){
                dtype tmp;
                if (!need_rand) {
                    in >> tmp;
                } else {
                    tmp = f(i, j); 
                }
                if (j % comm_size == rank) {
                    A(i, j - my_cols_offset) = tmp;
                    my_cols_offset += comm_size-1;

//                    if (j != glob_cols-1)
//                        tmp1[i] += tmp;
                }
            }
        }

//        MPI_Reduce(tmp1, sum_bvec, glob_rows, mpi_datatype, MPI_SUM, (glob_cols-1) % comm_size, MPI_COMM_WORLD);
//        if ((glob_cols-1) % comm_size == rank)
//            std::copy(sum_bvec, sum_bvec + glob_rows, &A(0, my_cols-1));

//         sleep(rand()/(float)RAND_MAX*5);
//         A.print();

        Matrix<dtype> Acopy = A;

        MPI_Barrier(MPI_COMM_WORLD);
        st_1 = MPI_Wtime();

        size_t size = glob_rows;
        dtype* x_vec = new dtype[size];

        my_cols_offset = rank;
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
                my_cols_offset += comm_size - 1;
            }
            MPI_Bcast(x_vec, size, mpi_datatype, root, MPI_COMM_WORLD);
            for (size_t col_i = ind/comm_size; col_i < my_cols; ++col_i) {
                dtype *y_vec = A.get_col(col_i);
                dtype scal = 2*scalar_prod(x_vec, y_vec, size);
                for (size_t j = 0; j < size; ++j) {
                    y_vec[j] -= scal*x_vec[j];
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        et_1 = MPI_Wtime();

        dtype* tmp_vec = new dtype[size+1];
        dtype* b_vec = new dtype[size];

        MPI_Barrier(MPI_COMM_WORLD);
        st_2 = MPI_Wtime();

        for (int i = 0; i < size; ++i) {
            x_vec[i] = 0;
            tmp_vec[i] = 0;
        }

        if (rank == (glob_cols - 1) % comm_size){
            dtype* ptr = A.get_col(A.n_cols() - 1);
            std::copy(ptr, ptr + size, b_vec);
        }
        MPI_Bcast(b_vec, size, mpi_datatype, (glob_cols - 1) % comm_size, MPI_COMM_WORLD);
//
//         sleep(rand()/(float)RAND_MAX*10);
//         A.print();

        if ((glob_cols-1) % comm_size == rank) {
            my_cols -= 1;
        }

        int offset = my_cols;
        for (int ind = size-1; ind >= 0; --ind){
            int root = ind % comm_size;

            dtype val = 0, val_rec = 0;
            for (int j = offset; j < my_cols; j += 1){

                val += A(ind, j)*x_vec[j*comm_size+rank];
//                std::cout << rank << "  "
//                          << ind << "  "
//                          << root << " | "
//                          << j << " / "
//                        << my_cols << " | "
//                        << x_vec[j*comm_size+rank]<< " | "
//                        << val << std::endl;
            }
            MPI_Reduce(&val, &val_rec, 1, mpi_datatype, MPI_SUM, root, MPI_COMM_WORLD);

            if (rank == root) {
                offset -= 1;
                x_vec[ind] = (b_vec[ind] - val_rec) / A(ind, ind/comm_size);
//                std::cout << "  " << rank << "  "
//                        << ind << "  "
//                        << b_vec[ind] << "  "
//                        << A(ind, ind/comm_size) << "  "
//                        << val_rec << "  "
//                          << x_vec[ind] << std::endl;
            } else {
                x_vec[ind] = 0;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        et_2 = MPI_Wtime();

        MPI_Allreduce(x_vec, tmp_vec, size, mpi_datatype, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < size; ++i) {
            x_vec[i] = tmp_vec[i];
        }

                delete[] b_vec;

        dtype* res_vec = new dtype[glob_rows];
        if (rank == (glob_cols-1) % comm_size) {
            b_vec = Acopy.get_col(Acopy.n_cols() - 1);
        }

        for (int i = 0; i < size; ++i){
            res_vec[i] = 0.0;
            tmp_vec[i] = 0.0;
            my_cols_offset = rank;
            for (int j = 0; j < my_cols; ++j){
                tmp_vec[i] += Acopy(i, j) * x_vec[my_cols_offset];
                my_cols_offset += comm_size;
            }
        }

        MPI_Reduce(tmp_vec, res_vec, size, mpi_datatype, MPI_SUM, (glob_cols-1) % comm_size, MPI_COMM_WORLD);
       
        if ((glob_cols-1) % comm_size == rank)
            for (size_t i = 0; i < size; ++i){
                //std::cout << res_vec[i] << " | " << b_vec[i] << std::endl;
                res_vec[i] -= b_vec[i];
            }

        dtype diff = 0.0;
        if ((glob_cols-1) % comm_size == rank)
            diff = norm(res_vec, size);
        
        long int elapsed_1, elapsed_1_tmp = ((et_1 - st_1) * 1000000);
        long int elapsed_2, elapsed_2_tmp = ((et_2 - st_2) * 1000000);

        
        MPI_Reduce(&elapsed_1_tmp, &elapsed_1, 1, mpi_datatype, MPI_MAX, (glob_cols-1) % comm_size, MPI_COMM_WORLD);
        MPI_Reduce(&elapsed_2_tmp, &elapsed_2, 1, mpi_datatype, MPI_MAX, (glob_cols-1) % comm_size, MPI_COMM_WORLD);

        if ((glob_cols-1) % comm_size == rank)
            for (size_t ind = 0; ind < size; ++ind) {
                //std::cout << "x_" << ind << ": " << x_vec[ind] << std::endl;
           }
        if ((glob_cols-1) % comm_size == rank){
            std::cout <<  "Mat_size " << size << " Comm_size " << comm_size << " Forward_Time_(microsec) " << elapsed_1 << "  Backward_Time_(microsec) " << elapsed_2 << " diff " << diff << std::endl;
        }
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
