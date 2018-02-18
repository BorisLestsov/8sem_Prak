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

        std::cout << "hi " << rank << std::endl;



    }
    catch (const std::string& e) {
        if (rank == 0)
            std::cerr << e << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch (const MPI::Exception& e) {
        if (rank == 0)
            std::cerr << e.Get_error_string() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Finalize();
    return 0;
}


