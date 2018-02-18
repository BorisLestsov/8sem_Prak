#include <iostream>
#include <fstream>
#include "Matrix.h"
#include <algorithm>
#include <vector>
#include <string>
#include "mpi.h"


using Matrix_ns::Matrix;

#define DEBUG

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

        Matrix<float> A(in, Matrix_ns::Normal, Matrix_ns::RowMaj);
        A.print();

        Matrix<float> U(4, 4, 0.0f);
        U(0,0) = 1;
        U(1,1) = 1;
        U(2,2) = 1;
        U(3,3) = 1;

        (A*U).print();

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


