#include <mpi.h>
#include <stdio.h>
#include <thread>
#include <string>
#include <chrono>
#include <ppltasks.h>
#include <windows.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

void writeVTK2(long timestep, bool* data, const char* prefix, int w, int h, int startx, int starty, int dim)
{
	char filename[2048];
	int x, y;

	int offsetX = 0;
	int offsetY = 0;
	float deltax = 1.0;
	long nxy = w * h * sizeof(float);

	snprintf(filename, sizeof(filename), "%s-%05ld%s", prefix, timestep, ".vti");
	FILE* fp;
	fopen_s(&fp, filename, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX,
	        offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
	fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

	for (y = starty; y < starty + dim; y++)
	{
		for (x = startx; x < startx + dim; x++)
		{
			float value = data[calcIndex(h, x,y)] ? 1.0f : 0.0f;
			fwrite((unsigned char*)&value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

void evolve(bool* currentfield, bool* newfield, int w, int h, int startx, int starty, int dim)
{
	int x, y;
	for (y = startx; y < startx + dim; y++)
	{
		for (x = starty; x < starty + dim; x++)
		{
			// DON'T NEST
			int nachbarnz = 0;
			for (int i = -1; i <= 1; ++i)
			{
				for (int j = -1; j <= 1; ++j)
				{
					if (x + i >= w || x + i <= 0 || y + j <= 0 || y + j >= h)
						continue;
					if (i == j && i == 0)
						continue;
					nachbarnz += currentfield[calcIndex(w, x + i, y + j)];
				}
			}

			bool wert = currentfield[calcIndex(w, x, y)];

			if (nachbarnz == 3 || (nachbarnz == 2 && wert))
			{
				newfield[calcIndex(w, x, y)] = true;
			}
			else
			{
				newfield[calcIndex(w, x, y)] = false;
			}
			//TODO FIXME impletent rules and assign new value

			//newfield[calcIndex(w, x, y)] = !currentfield[calcIndex(w, x, y)];
		}
	}
}

void doPrint(bool* field, FILE* file, int w, int startx, int starty, int dim)
{
	for (int y = startx; y < startx + dim; y++)
	{
		for (int x = starty; x < starty + dim; x++)
		{
			bool wert = field[calcIndex(w, x, y)];
			fwrite(&wert, 1, 1, file);
		}

		fwrite("\n", 1, 1, file);
	}
}

int _main(int argc, char* argv[], int* rank)
{
	auto now = std::chrono::high_resolution_clock::now();
	srand(std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count());

	int size, i;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const auto mySize = 15;
	const auto ppr = 2;
	const auto fieldSize = mySize * size;

	bool* field = static_cast<bool*>(calloc(fieldSize, sizeof(bool)));
	bool* newField = static_cast<bool*>(calloc(fieldSize, sizeof(bool)));
	for (int i = 0; i < fieldSize; ++i) field[i] = (rand() < RAND_MAX / 10) ? 1 : 0;

	int dim[2] = {ppr, ppr};
	int period[2] = {1, 1};
	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 1, &comm);

	int newRank;
	MPI_Comm_rank(comm, &newRank);
	*rank = newRank;

	int coords[2];
	MPI_Cart_coords(comm, newRank, 2, coords);
	coords[1] = ppr - coords[1] - 1;

	int left, right, up, down;
	MPI_Cart_shift(comm, 0, 1, &left, &right);
	MPI_Cart_shift(comm, 1, 1, &up, &down);


	int totalSize[2] = {ppr * mySize, ppr * mySize};

	int size_horizontal[2] = {1, mySize};
	int size_vertical[2] = {mySize, 1};

	// Inner datatypes
#pragma region inner
	int start_inner_left [2] = {coords[0] * mySize, coords[1] * mySize};
	start_inner_left[0] += totalSize[0];
	start_inner_left[0] %= totalSize[0];
	start_inner_left[1] += totalSize[1];
	start_inner_left[1] %= totalSize[1];

	MPI_Datatype inner_left;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_inner_left, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &inner_left);
	MPI_Type_commit(&inner_left);

	int start_inner_right [2] = {coords[0] * mySize + mySize - 1, coords[1] * mySize};
	start_inner_right[0] += totalSize[0];
	start_inner_right[0] %= totalSize[0];
	start_inner_right[1] += totalSize[1];
	start_inner_right[1] %= totalSize[1];
	MPI_Datatype inner_right;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_inner_right, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &inner_right);
	MPI_Type_commit(&inner_right);

	MPI_Datatype inner_up;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_inner_left, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &inner_up);
	MPI_Type_commit(&inner_up);

	int start_inner_down [2] = {coords[0] * mySize, coords[1] * mySize + mySize - 1};
	start_inner_down[0] += totalSize[0];
	start_inner_down[0] %= totalSize[0];
	start_inner_down[1] += totalSize[1];
	start_inner_down[1] %= totalSize[1];

	MPI_Datatype inner_down;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_inner_down, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &inner_down);
	MPI_Type_commit(&inner_down);
#pragma endregion

	int start_outer_left [2] = {coords[0] * mySize - 1, coords[1] * mySize};
	start_outer_left[0] += totalSize[0];
	start_outer_left[0] %= totalSize[0];
	start_outer_left[1] += totalSize[1];
	start_outer_left[1] %= totalSize[1];
	MPI_Datatype outer_left;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_inner_left, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &outer_left);
	MPI_Type_commit(&outer_left);

	int start_outer_right [2] = {coords[0] * mySize + mySize + 1, coords[1] * mySize};
	start_outer_right[0] += totalSize[0];
	start_outer_right[0] %= totalSize[0];
	start_outer_right[1] += totalSize[1];
	start_outer_right[1] %= totalSize[1];
	MPI_Datatype outer_right;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_outer_right, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &outer_right);
	MPI_Type_commit(&outer_right);

	int start_outer_up [2] = {coords[0] * mySize, coords[1] * mySize - 1};
	start_outer_up[0] += totalSize[0];
	start_outer_up[0] %= totalSize[0];
	start_outer_up[1] += totalSize[1];
	start_outer_up[1] %= totalSize[1];
	MPI_Datatype outer_up;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_outer_up, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &outer_up);
	MPI_Type_commit(&outer_up);

	int start_outer_down [2] = {coords[0] * mySize, coords[1] * mySize + mySize + 1};
	start_outer_down[0] += totalSize[0];
	start_outer_down[0] %= totalSize[0];
	start_outer_down[1] += totalSize[1];
	start_outer_down[1] %= totalSize[1];
	MPI_Datatype outer_down;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_outer_down, MPI_ORDER_C, MPI_UNSIGNED_CHAR,
	                         &outer_down);
	MPI_Type_commit(&outer_down);

	const auto print = [&](bool* field, int t, int startx, int starty, int dim)
	{
		std::string name = "gol-" + std::to_string(newRank) + "-" + std::to_string(t) + ".rip";
		FILE* f;
		fopen_s(&f, name.data(), "wb");

		doPrint(field, f, ppr * mySize, startx, starty, dim);

		fclose(f);
	};

	printf("[%d] before-loop\n", newRank);
	fflush(stdout);

	int t = 0;
	while (true)
	{
		printf("[%d] begin-loop\n", newRank);
		fflush(stdout);

		MPI_Request req[16];
		MPI_Request* req_ptr = req;

		MPI_Irecv(field, 1, outer_right, right, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_left, left, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_up, up, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_down, down, 1, comm, req_ptr++);

		MPI_Isend(field, 1, inner_left, left, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_right, right, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_up, up, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_down, down, 1, comm, req_ptr++);

		MPI_Status _status[16];
		MPI_Waitall(8, req, _status);

		if (t == 3)
		{
			printf("[%d] exitting\n", newRank);
			fflush(stdout);
			break;
		}

		printf("[%d] before-print\n", newRank);
		fflush(stdout);


		//print(field, t, start_inner_left[0], start_inner_left[1], mySize);

		std::string name = "gol-" + std::to_string(newRank) + "-";
		writeVTK2(t, field, name.data(), ppr * mySize, ppr * mySize, start_inner_left[0], start_inner_left[1], mySize);

		printf("[%d] after-print\n", newRank);
		fflush(stdout);
		evolve(field, newField, ppr * mySize, ppr * mySize, start_inner_left[0], start_inner_left[1], mySize);
		printf("[%d] after-evolve\n", newRank);
		fflush(stdout);
		std::swap(field, newField);
		++t;
	}

	MPI_Finalize();

	return 0;
}

int rank = -1;

LONG __stdcall filter(_EXCEPTION_POINTERS* p)
{
	printf("[%d] Crashed!\n", rank);
	fflush(stdout);
	return EXCEPTION_CONTINUE_SEARCH;
}

int main(int argc, char* argv[])
{
	SetUnhandledExceptionFilter(filter);
	AddVectoredExceptionHandler(TRUE, filter);
	_main(argc, argv, &rank);
}
