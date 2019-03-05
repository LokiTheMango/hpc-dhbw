#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <time.h>
#else
#include <endian.h>
#include <unistd.h>
#include <sys/time.h>
#endif
#include <stdio.h>
//#include <cstdlib.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <fstream>
#include <thread>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

long TimeSteps = 100;

void writeVTK2(long timestep, double* data, const char* prefix, int w, int h)
{
	char filename[2048];
	int x, y;

	int offsetX = 0;
	int offsetY = 0;
	float deltax = 1.0;
	long nxy = w * h * sizeof(float);

	snprintf(filename, sizeof(filename), "%s-%05ld%s", prefix, timestep, ".vti");
	FILE* fp = fopen(filename, "w");

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

	for (y = 0; y < h; y++)
	{
		for (x = 0; x < w; x++)
		{
			float value = data[calcIndex(h, x,y)];
			fwrite((unsigned char*)&value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}


void show(double* currentfield, int w, int h)
{
	printf("\033[H");
	int x, y;
	for (y = 0; y < h; y++)
	{
		for (x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
		//printf("\033[E");
		printf("\n");
	}
	fflush(stdout);
}


void evolve(double* currentfield, double* newfield, int w, int h)
{
	int x, y;
	for (y = 0; y < h; y++)
	{
		for (x = 0; x < w; x++)
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

			double wert = currentfield[calcIndex(w, x, y)];

			if (nachbarnz == 3 || (nachbarnz == 2 && wert))
			{
				newfield[calcIndex(w, x, y)] = 1.0;
			}
			else
			{
				newfield[calcIndex(w, x, y)] = 0;
			}
			//TODO FIXME impletent rules and assign new value

			//newfield[calcIndex(w, x, y)] = !currentfield[calcIndex(w, x, y)];
		}
	}
}

void filling(double* currentfield, int w, int h)
{
	int i;
	for (i = 0; i < h * w; i++)
	{
		currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
	}
}

void game(int w, int h)
{
	int rank, size;
	double* currentfield = (double*)calloc(w * h, sizeof(double));
	double* newfield = (double*)calloc(w * h, sizeof(double));

	//printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

	filling(currentfield, w, h);
	long t;
	for (t = 0;; t++)
	{
		show(currentfield, w, h);
		evolve(currentfield, newfield, w, h);

		printf("%ld timestep\n", t);
		//writeVTK2(t, currentfield, "gol", w, h);

		using namespace std::literals;
		std::this_thread::sleep_for(200ms);

		//SWAP
		double* temp = currentfield;
		currentfield = newfield;
		newfield = temp;

		if(!std::memcmp(currentfield, newfield, w * h * sizeof(double))) break;
	}

	free(currentfield);
	free(newfield);
}

int ___main(int c, char** v)
{
	int w = 0, h = 0;
	if (c > 1) w = atoi(v[1]); ///< read width
	if (c > 2) h = atoi(v[2]); ///< read height
	if (w <= 0) w = 30; ///< default width
	if (h <= 0) h = 30; ///< default height
	game(w, h);

	return 0;
}
