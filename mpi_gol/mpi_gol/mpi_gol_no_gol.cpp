#include <mpi.h>
#include <stdio.h>
#include <thread>
#include <string>
#include <chrono>

int ___main(int c, char** v);

int main(int argc, char* argv[])
{
	auto now = std::chrono::high_resolution_clock::now();
	srand(std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count());

	int rank, size, i;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const auto mySize = 1;
	const auto fieldSize = mySize * size;

	int* field = static_cast<int*>(calloc(fieldSize, sizeof(int)));
	for(int i = 0; i < fieldSize; ++i) field[i] = i;

	int dim = size;
	int period = 1;
	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, 1, &dim, &period, 1, &comm);

	int newRank;
	MPI_Comm_rank(comm, &newRank);

	int src, dst;
	MPI_Cart_shift(comm, 0, 1, &src, &dst);


	int lsize = 1;
	int start = newRank * mySize - 1 + fieldSize;
	start %= fieldSize;
	MPI_Datatype left;
	MPI_Type_create_subarray(1, &fieldSize, &lsize, &start, MPI_ORDER_C, MPI_INT, &left);
	MPI_Type_commit(&left);

	int start2 = newRank * mySize + 1 + fieldSize;
	start2 %= fieldSize;
	MPI_Datatype right;
	MPI_Type_create_subarray(1, &fieldSize, &lsize, &start2, MPI_ORDER_C, MPI_INT, &right);
	MPI_Type_commit(&right);

	int start3 = newRank * mySize + fieldSize;
	start3 %= fieldSize;
	MPI_Datatype inner;
	MPI_Type_create_subarray(1, &fieldSize, &lsize, &start3, MPI_ORDER_C, MPI_INT, &inner);
	MPI_Type_commit(&inner);

	field[newRank * mySize] = rand();

	MPI_Request req[4];
	MPI_Irecv(field, 1, right, dst, 1, comm, req);
	MPI_Irecv(field, 1, left, src, 2, comm, req + 1);
	MPI_Isend(field, 1, inner, dst, 2, comm, req + 2);
	MPI_Isend(field, 1, inner, src, 1, comm, req + 3);

	MPI_Status _status[4];
	MPI_Waitall(4, req, _status);

	std::string _string;
	_string += "Rank:\t" + std::to_string(newRank) + " - " + std::to_string(src) + " " + std::to_string(dst) + "\n";
	_string += "My: " + std::to_string(field[start3]) + "\n";
	_string += "Left:\t" + std::to_string(field[start]) + "\n";
	_string += "Right:\t" + std::to_string(field[start2]) + "\n\n";

	printf("%s", _string.data());

	MPI_Finalize();

	free(field);
	return 0;
}
