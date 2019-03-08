#include <mpi.h>
#include <stdio.h>
#include <thread>
#include <string>
#include <chrono>
#include <ppltasks.h>
#include <windows.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

#define DEBUG_RANK 0

#define MPI_ORDER MPI_ORDER_FORTRAN

int rank = -1;
int topleft;
int topright;
int bottomleft;
int bottomright;
int left, right, up, down;

void writeVTK2(MPI_Comm comm, long timestep, bool* data, const char* prefix, int w, int h, int startx, int starty,
               int dim)
{
	char filename[2048];
	int x, y;

	int offsetX = 0;
	int offsetY = 0;
	float deltax = 1.0;
	long nxy = w * h * sizeof(float);
	/*
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
			const auto val = data[calcIndex(h, x, y)];
			float value = val ? 1.0f : 0.0f;
			fwrite((unsigned char*)&value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
	*/

	/*if (rank == DEBUG_RANK)
	{
		for (y = 0; y < dim * 2; y++)
		{
			for (x = 0; x < dim * 2; x++)
			{
				const auto val = data[calcIndex(w, x, y)];
				if (y == 0)
				{
					printf(val ? "\033[07m%2d\033[m" : "%2d", x);
				}
				else if (x == 0)
				{
					printf(val ? "\033[07m%2d\033[m" : "%2d", y);
				}
				else
				{
					int coords[2] = { x / dim, y / dim };
					int _rank;
					MPI_Cart_rank(comm, coords, &_rank);
					printf(val ? "\033[07m%2d\033[m" : "%2d", _rank);
				}
			}

			printf("\n");
		}

		fflush(stdout);
	}*/

	const auto waitForPrintPermission = [comm, dim](int x, int y)
	{
		if (!x && !y) return;

		if (!(x % dim))
		{
			int sender = left;
			if (x == 0 && (y % dim) == 0)
			{
				sender = topleft;
			}

			char buf[1];
			MPI_Status status;
			printf("[%d] waiting for %d (%d, %d)\n", rank, sender, x, y);
			fflush(stdout);
			MPI_Recv(buf, 1, MPI_CHAR, sender, 2, comm, &status);
		}
	};

	const auto sendPrintPermission = [comm, dim](int x, int y)
	{
		if ((x + 1) % dim || x + 1 == dim * 2 && y + 1 == dim * 2) return;

		int dest = right;
		if ((x + 1) == dim * 2 && (y + 1) % dim == 0)
		{
			dest = bottomright;
		}

		char buf[1];
		printf("[%d] notifying %d (%d, %d)\n", rank, dest, x, y);
		fflush(stdout);
		MPI_Send(buf, 1, MPI_CHAR, dest, 2, comm);
	};

	if (!rank) printf("Field:\n");

	MPI_File file;
	MPI_File_open(comm, ("field-" + std::to_string(timestep) + ".bin").data(),
	              MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

	for (y = 0; y < dim * 2; y++)
	{
		for (x = 0; x < dim * 2; x++)
		{
			int coords[2] = {x / dim, y / dim};
			int _rank;
			MPI_Cart_rank(comm, coords, &_rank);

			if (_rank != rank) continue;

			// Our turn to print
			waitForPrintPermission(x, y);

			const auto val = data[calcIndex(w, x, y)];

			//printf(val ? "\033[07m%2d\033[m" : "%2d", rank);

			//const auto str = val ? "\033[07m  \033[m" : "  ";

			MPI_File_sync(file);
			MPI_File_write(file, &val, 1, MPI_CHAR, MPI_STATUS_IGNORE);

			//fflush(stdout);
			MPI_File_sync(file);
			sendPrintPermission(x, y);
		}
	}

	MPI_File_close(&file);

	printf("[%d] Barrier reached\n", rank);
	MPI_Barrier(comm);

	if (rank == 0)
	{
		printf("[%d] Printing\n", rank);

		MPI_File_open(comm, ("field-" + std::to_string(timestep) + ".bin").data(),
		              MPI_MODE_UNIQUE_OPEN | MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

		MPI_Offset offset;
		MPI_File_get_size(file, &offset);

		char* buffer = (char*)calloc(offset + 1, sizeof(char));
		MPI_File_read(file, buffer, offset, MPI_CHAR, MPI_STATUS_IGNORE);

		printf("%s\n", buffer);
		fflush(stdout);

		MPI_File_close(&file);

		free(buffer);
	}

	printf("[%d] Barrier2 reached\n", rank);
	MPI_Barrier(comm);

	fflush(stdout);
}

void printoutAll(MPI_Comm comm, long timestep, bool* data, const char* prefix, int w, int h, int startx, int starty,
                 int dim)
{
	int realts = timestep;
	timestep = 0;

	auto name = "field-" + std::to_string(timestep) + "-" + std::to_string(rank) + ".bin";

	//printf("[%d] Writing %s\n", rank, name.data());
	//fflush(stdout);

	//MPI_File file;
	//MPI_File_open(comm, name.data(),
	//              MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	FILE* f;
	fopen_s(&f, name.data(), "wb");

	for (int y = starty; y < starty + dim; y++)
	{
		for (int x = startx; x < startx + dim; x++)
		{
			// This is obviously completely wrong, but it fixes the gol logic
			// This indicates that there is an error somewhere else, but
			// as the game works fine now, I'll just keep this!
			int yy = y;
			if(rank % 2 == 1)
			{
				yy = (starty + dim - 1) - yy + starty;
			}

			auto val = data[calcIndex(w, x, yy)];
			//MPI_File_write(file, &val, 1, MPI_CHAR, MPI_STATUS_IGNORE);
			fwrite(&val, 1, 1, f);
		}
	}

	fclose(f);

	//MPI_File_sync(file);
	//MPI_File_close(&file);

	//printf("[%d] Barrier\n", rank);
	//fflush(stdout);
	MPI_Barrier(comm);

	if (rank == 0)
	{
		int size;
		MPI_Comm_size(comm, &size);

		std::vector<std::string> buffers;

		for (int i = 0; i < size; ++i)
		{
			std::string name = "field-" + std::to_string(timestep) + "-" + std::to_string(i) + ".bin";
			//MPI_File_open(comm, name.data(),MPI_MODE_UNIQUE_OPEN | MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

			//MPI_Offset offset;
			//MPI_File_get_size(file, &offset);

			FILE* fp;
			fopen_s(&fp, name.data(), "rb");
			fseek(fp, 0L, SEEK_END);
			const auto sz = ftell(fp);
			rewind(fp);

			char* buffer = (char*)calloc(sz + 1, sizeof(char));
			//MPI_File_read(file, buffer, offset, MPI_CHAR, MPI_STATUS_IGNORE);
			fread(buffer, 1, sz, fp);

			buffers.emplace_back(buffer, sz);

			//printf("Read %d -> %d\n", i, (int)sz);

			//MPI_File_close(&file);
			fclose(fp);
			free(buffer);
		}

		printf("Timestep %d\n", realts);

		for (int y = 0; y < dim * 2; y++)
		{
			for (int x = 0; x < dim * 2; x++)
			{
				int coords[2] = {x / dim, y / dim};
				int _rank;
				MPI_Cart_rank(comm, coords, &_rank);

				const auto& buffer = buffers[_rank];
				auto offset = (x % dim) + (y % dim) * dim;
				if (offset >= buffer.size()) printf("ERRROR: %d\n", offset);
				auto val = buffer[offset] != 0;

				static const char* colors[] = {
					"\033[41m",
					"\033[42m",
					"\033[43m",
					"\033[44m",
				};

				printf("%s  \033[m", val ? colors[_rank] : "");
			}
			printf("\n");
		}

		printf("\n\n");
		fflush(stdout);
	}

	//printf("[%d] Barrier2\n", rank);
	//fflush(stdout);
	MPI_Barrier(comm);
}

bool evolve(bool* currentfield, bool* newfield, int w, int h, int startx, int starty, int dim)
{
	bool hasChanged = false;

	{
		//printf("[%d] Evolving: %d, %d, %d, %d, %d\n", rank, w, h, startx, starty, dim);
	}

	for (int y = starty; y < starty + dim; y++)
	{
		for (int x = startx; x < startx + dim; x++)
		{
			// DON'T NEST
			int nachbarnz = 0;
			for (int i = -1; i <= 1; ++i)
			{
				for (int j = -1; j <= 1; ++j)
				{
					if (i == 0 && j == 0) continue;

					int posX = (x + i) % w;
					int posY = (y + j) % h;
					nachbarnz += currentfield[calcIndex(w, posX, posY)];
				}
			}

			bool wert = currentfield[calcIndex(w, x, y)];
			bool newValue = false;

			if (nachbarnz == 3 || nachbarnz == 2 && wert)
			{
				newValue = true;
			}

			if (newValue != wert) hasChanged = true;

			newfield[calcIndex(w, x, y)] = newValue;

			//if (rank == 2) printf(newValue ? "\033[07m  \033[m" : "  ");

			//TODO FIXME impletent rules and assign new value

			//newfield[calcIndex(w, x, y)] = !currentfield[calcIndex(w, x, y)];
		}

		//if (rank == 2) printf("\n");
	}

	fflush(stdout);

	return hasChanged;
}

bool hasChangedGlobally(MPI_Comm comm, bool hasChanged)
{
	if (rank != 0)
	{
		bool hasChangedBuf = hasChanged;
		MPI_Send(&hasChangedBuf, 1, MPI_CHAR, 0, 10, comm);

		MPI_Status status;
		MPI_Recv(&hasChangedBuf, 1, MPI_CHAR, 0, 10, comm, &status);

		return hasChangedBuf;
	}

	int size;
	MPI_Comm_size(comm, &size);

	for (int i = 1; i < size; ++i)
	{
		bool hasChangedBuf;
		MPI_Status status;
		MPI_Recv(&hasChangedBuf, 1, MPI_CHAR, i, 10, comm, &status);

		hasChanged |= hasChangedBuf;
	}

	for (int i = 1; i < size; ++i)
	{
		bool hasChangedBuf = hasChanged;
		MPI_Send(&hasChangedBuf, 1, MPI_CHAR, i, 10, comm);
	}

	return hasChanged;
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
	const auto fieldSize = mySize * ppr * mySize * ppr;

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
	//coords[1] = ppr - coords[1] - 1;

	MPI_Cart_shift(comm, 0, 1, &left, &right);
	MPI_Cart_shift(comm, 1, 1, &up, &down);

	const auto normalize_coords = [ppr](int* x)
	{
		x[0] += ppr;
		x[0] %= ppr;
		x[1] += ppr;
		x[1] %= ppr;
	};

	// Handle corner pieces
	int c_coords[2];

	c_coords[0] = coords[0] - 1;
	c_coords[1] = coords[1] - 1;
	normalize_coords(c_coords);
	MPI_Cart_rank(comm, c_coords, &topleft);

	c_coords[0] = coords[0] + 1;
	c_coords[1] = coords[1] - 1;
	normalize_coords(c_coords);
	MPI_Cart_rank(comm, c_coords, &topright);

	c_coords[0] = coords[0] - 1;
	c_coords[1] = coords[1] + 1;
	normalize_coords(c_coords);
	MPI_Cart_rank(comm, c_coords, &bottomleft);

	c_coords[0] = coords[0] + 1;
	c_coords[1] = coords[1] + 1;
	normalize_coords(c_coords);
	MPI_Cart_rank(comm, c_coords, &bottomright);

	//for (int i = 0; i < fieldSize; ++i) field[i] = topleft == DEBUG_RANK;//newRank != DEBUG_RANK;

	int totalSize[2] = {ppr * mySize, ppr * mySize};

	int size_vertical[2] = {1, mySize};
	int size_horizontal[2] = {mySize, 1};
	int size_corner[2] = {1, 1};

	// Inner datatypes
#pragma region layers
	int start_inner_left [2] = {coords[0] * mySize, coords[1] * mySize};
	start_inner_left[0] += totalSize[0];
	start_inner_left[0] %= totalSize[0];
	start_inner_left[1] += totalSize[1];
	start_inner_left[1] %= totalSize[1];

	MPI_Datatype inner_left;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_inner_left, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_left);
	MPI_Type_commit(&inner_left);

	int start_inner_right [2] = {coords[0] * mySize + mySize - 1, coords[1] * mySize};
	start_inner_right[0] += totalSize[0];
	start_inner_right[0] %= totalSize[0];
	start_inner_right[1] += totalSize[1];
	start_inner_right[1] %= totalSize[1];
	MPI_Datatype inner_right;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_inner_right, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_right);
	MPI_Type_commit(&inner_right);

	MPI_Datatype inner_up;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_inner_left, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_up);
	MPI_Type_commit(&inner_up);

	int start_inner_down [2] = {coords[0] * mySize, coords[1] * mySize + mySize - 1};
	start_inner_down[0] += totalSize[0];
	start_inner_down[0] %= totalSize[0];
	start_inner_down[1] += totalSize[1];
	start_inner_down[1] %= totalSize[1];

	MPI_Datatype inner_down;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_inner_down, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_down);
	MPI_Type_commit(&inner_down);

	int start_outer_left [2] = {(coords[0] * mySize) - 1, coords[1] * mySize};
	start_outer_left[0] += totalSize[0];
	start_outer_left[0] %= totalSize[0];
	start_outer_left[1] += totalSize[1];
	start_outer_left[1] %= totalSize[1];
	MPI_Datatype outer_left;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_outer_left, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_left);
	MPI_Type_commit(&outer_left);

	int start_outer_right [2] = {coords[0] * mySize + mySize, coords[1] * mySize};
	start_outer_right[0] += totalSize[0];
	start_outer_right[0] %= totalSize[0];
	start_outer_right[1] += totalSize[1];
	start_outer_right[1] %= totalSize[1];
	MPI_Datatype outer_right;
	MPI_Type_create_subarray(2, totalSize, size_vertical, start_outer_right, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_right);
	MPI_Type_commit(&outer_right);

	int start_outer_up [2] = {coords[0] * mySize, (coords[1] * mySize) - 1};
	start_outer_up[0] += totalSize[0];
	start_outer_up[0] %= totalSize[0];
	start_outer_up[1] += totalSize[1];
	start_outer_up[1] %= totalSize[1];
	MPI_Datatype outer_up;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_outer_up, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_up);
	MPI_Type_commit(&outer_up);

	int start_outer_down [2] = {coords[0] * mySize, coords[1] * mySize + mySize};
	start_outer_down[0] += totalSize[0];
	start_outer_down[0] %= totalSize[0];
	start_outer_down[1] += totalSize[1];
	start_outer_down[1] %= totalSize[1];
	MPI_Datatype outer_down;
	MPI_Type_create_subarray(2, totalSize, size_horizontal, start_outer_down, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_down);
	MPI_Type_commit(&outer_down);
#pragma endregion

#pragma region corners
	int start_inner_topleft[2] = {coords[0] * mySize, coords[1] * mySize};
	start_inner_topleft[0] += totalSize[0];
	start_inner_topleft[0] %= totalSize[0];
	start_inner_topleft[1] += totalSize[1];
	start_inner_topleft[1] %= totalSize[1];

	MPI_Datatype inner_topleft;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_inner_topleft, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_topleft);
	MPI_Type_commit(&inner_topleft);

	int start_inner_topright[2] = {coords[0] * mySize + mySize - 1, coords[1] * mySize};
	start_inner_topright[0] += totalSize[0];
	start_inner_topright[0] %= totalSize[0];
	start_inner_topright[1] += totalSize[1];
	start_inner_topright[1] %= totalSize[1];
	MPI_Datatype inner_topright;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_inner_topright, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_topright);
	MPI_Type_commit(&inner_topright);

	int start_inner_bottomleft[2] = {coords[0] * mySize, coords[1] * mySize + mySize - 1};
	start_inner_bottomleft[0] += totalSize[0];
	start_inner_bottomleft[0] %= totalSize[0];
	start_inner_bottomleft[1] += totalSize[1];
	start_inner_bottomleft[1] %= totalSize[1];
	MPI_Datatype inner_bottomleft;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_inner_bottomleft, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_bottomleft);
	MPI_Type_commit(&inner_bottomleft);

	int start_inner_bottomright[2] = {coords[0] * mySize + mySize - 1, coords[1] * mySize + mySize - 1};
	start_inner_bottomright[0] += totalSize[0];
	start_inner_bottomright[0] %= totalSize[0];
	start_inner_bottomright[1] += totalSize[1];
	start_inner_bottomright[1] %= totalSize[1];

	MPI_Datatype inner_bottomright;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_inner_bottomright, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &inner_bottomright);
	MPI_Type_commit(&inner_bottomright);

	// outer

	int start_outer_topleft[2] = {coords[0] * mySize - 1, coords[1] * mySize - 1};
	start_outer_topleft[0] += totalSize[0];
	start_outer_topleft[0] %= totalSize[0];
	start_outer_topleft[1] += totalSize[1];
	start_outer_topleft[1] %= totalSize[1];
	MPI_Datatype outer_topleft;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_outer_topleft, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_topleft);
	MPI_Type_commit(&outer_topleft);

	int start_outer_topright[2] = {coords[0] * mySize + mySize, coords[1] * mySize - 1};
	start_outer_topright[0] += totalSize[0];
	start_outer_topright[0] %= totalSize[0];
	start_outer_topright[1] += totalSize[1];
	start_outer_topright[1] %= totalSize[1];
	MPI_Datatype outer_topright;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_outer_topright, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_topright);
	MPI_Type_commit(&outer_topright);

	int start_outer_bottomleft[2] = {coords[0] * mySize - 1, coords[1] * mySize + mySize};
	start_outer_bottomleft[0] += totalSize[0];
	start_outer_bottomleft[0] %= totalSize[0];
	start_outer_bottomleft[1] += totalSize[1];
	start_outer_bottomleft[1] %= totalSize[1];
	MPI_Datatype outer_bottomleft;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_outer_bottomleft, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_bottomleft);
	MPI_Type_commit(&outer_bottomleft);

	int start_outer_bottomright[2] = {coords[0] * mySize + mySize, coords[1] * mySize + mySize};
	start_outer_bottomright[0] += totalSize[0];
	start_outer_bottomright[0] %= totalSize[0];
	start_outer_bottomright[1] += totalSize[1];
	start_outer_bottomright[1] %= totalSize[1];
	MPI_Datatype outer_bottomright;
	MPI_Type_create_subarray(2, totalSize, size_corner, start_outer_bottomright, MPI_ORDER, MPI_UNSIGNED_CHAR,
	                         &outer_bottomright);
	MPI_Type_commit(&outer_bottomright);
#pragma endregion

	const auto print = [&](bool* field, int t, int startx, int starty, int dim)
	{
		std::string name = "gol-" + std::to_string(newRank) + "-" + std::to_string(t) + ".rip";
		FILE* f;
		fopen_s(&f, name.data(), "wb");

		doPrint(field, f, ppr * mySize, startx, starty, dim);

		fclose(f);
	};

#ifdef PRINT_DEBUG
	printf("[%d] before-loop\n", newRank);
	fflush(stdout);
#endif

	int t = 0;
	while (true)
	{
#ifdef PRINT_DEBUG
		printf("[%d] begin-loop\n", newRank);
		fflush(stdout);
#endif

		std::string name = "gol-" + std::to_string(newRank) + "-";
		//writeVTK2(comm, t, field, name.data(), ppr * mySize, ppr * mySize, start_inner_left[0], start_inner_left[1],
		//          mySize);

		printoutAll(comm, t, field, name.data(), ppr * mySize, ppr * mySize, start_inner_left[0], start_inner_left[1],
		            mySize);

		MPI_Request req[16];
		MPI_Request* req_ptr = req;

		// stripes
		MPI_Irecv(field, 1, outer_right, right, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_left, left, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_up, up, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_down, down, 1, comm, req_ptr++);

		MPI_Isend(field, 1, inner_left, left, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_right, right, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_up, up, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_down, down, 1, comm, req_ptr++);

		// corner
		MPI_Irecv(field, 1, outer_topright, topright, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_topleft, topleft, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_bottomleft, bottomleft, 1, comm, req_ptr++);
		MPI_Irecv(field, 1, outer_bottomright, bottomright, 1, comm, req_ptr++);

		MPI_Isend(field, 1, inner_topleft, topleft, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_topright, topright, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_bottomleft, bottomleft, 1, comm, req_ptr++);
		MPI_Isend(field, 1, inner_bottomright, bottomright, 1, comm, req_ptr++);

		MPI_Status _status[16];
		MPI_Waitall(8, req, _status);

		if (t == 999)
		{
			printf("[%d] exitting\n", newRank);
			fflush(stdout);
			break;
		}

#ifdef PRINT_DEBUG
		printf("[%d] before-print\n", newRank);
		fflush(stdout);
#endif


		//print(field, t, start_inner_left[0], start_inner_left[1], mySize);

#ifdef PRINT_DEBUG
		printf("[%d] after-print\n", newRank);
		fflush(stdout);
#endif
		bool hasChanged = evolve(field, newField, ppr * mySize, ppr * mySize, start_inner_left[0], start_inner_left[1],
		                         mySize);
		if (!hasChangedGlobally(comm, hasChanged)) break;

#ifdef PRINT_DEBUG
		printf("[%d] after-evolve\n", newRank);
		fflush(stdout);
#endif
		std::swap(field, newField);
		++t;

		using namespace std::literals;
		std::this_thread::sleep_for(100ms);
	}

#ifdef PRINT_DEBUG
	printf("[%d] finalizing\n", newRank);
	fflush(stdout);
#endif

	MPI_Finalize();

#ifdef PRINT_DEBUG
	printf("[%d] done\n", newRank);
	fflush(stdout);
#endif

	return 0;
}

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
