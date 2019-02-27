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
#include <stdlib.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <fstream>

#define min(x, y) (x > y ? y : x)

#define calcIndex(width, x, y) ((y) * (width) + (x))

long TimeSteps = 100;

void writeVTK2(long timestep, bool* data, const char* prefix, int w, int h, int offsetX, int offsetY)
{
	char filename[2048];
	int x, y;

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
	fwrite((unsigned char *)&nxy, sizeof(long), 1, fp);

	for (y = 0; y < h; y++)
	{
		for (x = 0; x < w; x++)
		{
			float value = data[calcIndex(h, x, y)] ? 1.0f : 0.0f;
			fwrite((unsigned char *)&value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

void show(bool** currentfield, int w, int h, int num, int xperrow)
{
	printf("\033[H");

	int x, y;
	for (y = 0; y < h * xperrow; y++)
	{
		for (x = 0; x < w * xperrow; x++)
		{
			int localyIndex = (y % h);
			int localxIndex = (x % w);

			int yIndex = (y / h);
			int totalIndex = (yIndex * xperrow) + (x / w);
			const auto field = currentfield[totalIndex];

			bool val = field[calcIndex(w + 2, localxIndex + 1, localyIndex + 1)];

			if (val)
			{
				printf("");
			}

			printf(val ? "\033[07m  \033[m" : "  ");
		}
		printf("\n");
	}
	fflush(stdout);
}

void evolve(bool* currentfield, bool* newfield, int w, int h)
{
	int x, y;
	for (y = 1; y < h + 1; y++)
	{
		for (x = 1; x < w + 1; x++)
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
					nachbarnz += currentfield[calcIndex(w + 2, x + i, y + j)];
				}
			}

			bool wert = currentfield[calcIndex(w + 2, x, y)];

			if (nachbarnz == 3 || (nachbarnz == 2 && wert))
			{
				newfield[calcIndex(w + 2, x, y)] = 1;
			}
			else
			{
				newfield[calcIndex(w + 2, x, y)] = 0;
			}
			//TODO FIXME impletent rules and assign new value

			//newfield[calcIndex(w, x, y)] = !currentfield[calcIndex(w, x, y)];
		}
	}
}

void filling(bool* currentfield, int w, int h)
{
	int i;
	for (i = 0; i < h * w; i++)
	{
		currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
	}
	return;

	std::ifstream im("im.txt");
	im.seekg(0, std::ios::end);
	auto len = im.tellg();
	im.seekg(0, std::ios::beg);


	std::string field;
	field.resize(len);
	im.read((char*)field.data(), len);
	im.close();

	int pos = 0;
	for (int i = 0; i < h * w; ++i)
	{
		char x = field[pos++ % field.size()];
		if (x == '0')
		{
			currentfield[i] = 0.0;
		}
		else if (x == '1')
		{
			currentfield[i] = 1.0;
		}
		else
		{
			--i;
		}
	}
}

void synchronizeBorders(bool** fields, int num, int xPerrow, int w, int h)
{
	int rowSize = w + 2;
	int height = h + 2;

	for (int i = 0; i < num; ++i)
	{
		bool* field = fields[i];

		int leftIndex = i - 1;
		int rightIndex = i + 1;
		int upIndex = i - xPerrow;
		int downIndex = i + xPerrow;

		if (leftIndex >= 0)
		{
			bool* leftField = fields[leftIndex];

			for (int y = 0; y < height; ++y)
			{
				field[y * rowSize] = leftField[(y + 1) * rowSize - 1];
			}
		}

		if (upIndex >= 0)
		{
			bool* upField = fields[upIndex];

			for (int x = 0; x < rowSize; ++x)
			{
				field[x] = upField[x + (height - 1) * rowSize];
			}
		}

		if (rightIndex < num)
		{
			bool* rightField = fields[rightIndex];

			for (int y = 0; y < height; ++y)
			{
				field[(y + 1) * rowSize - 1] = rightField[y * rowSize];
			}
		}

		if (downIndex < num)
		{
			bool* downField = fields[downIndex];

			for (int x = 0; x < rowSize; ++x)
			{
				field[x + (height - 1) * rowSize] = downField[x];
			}
		}
	}
}

void game(int w, int h)
{
	bool activeField = false;

	//printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

	int numThreads = 4;
	omp_set_num_threads(numThreads);

	int threadnums = numThreads; //omp_get_num_threads();

	static const int fildDim = 15;

	static const int xPreRow = 2;
	static const int realRowSize = fildDim;
	static const int rowsize = realRowSize + 2;

	static const int realfieldsize = realRowSize * realRowSize;
	static const int fieldsize = rowsize * rowsize;

	bool** fields = (bool* *)calloc(threadnums, sizeof(bool*));
	bool** otherfields = (bool* *)calloc(threadnums, sizeof(bool*));

	for (int i = 0; i < threadnums; ++i)
	{
		fields[i] = (bool*)calloc(fieldsize, sizeof(bool));
		otherfields[i] = (bool*)calloc(fieldsize, sizeof(bool));
		filling(fields[i], rowsize, rowsize);
		filling(otherfields[i], rowsize, rowsize);
	}

#pragma omp parallel shared(threadnums, fields, otherfields, activeField, xPreRow, fildDim)
	{
		long t;
		for (t = 0; ; t++)
		{
#pragma omp master
			{
				synchronizeBorders((activeField ? fields : otherfields), threadnums, xPreRow, fildDim, fildDim);
				show((activeField ? fields : otherfields), fildDim, fildDim, threadnums, xPreRow);
			}

#pragma omp barrier

			{
				int threadnum = omp_get_thread_num();

				bool* mycurrentfield = (activeField ? fields : otherfields)[threadnum];
				bool* mynewfield = (!activeField ? fields : otherfields)[threadnum];

				evolve(mycurrentfield, mynewfield, fildDim, fildDim);

				printf("%ld timestep\n", t);

				std::string name("gol-" + std::to_string(omp_get_thread_num()));
				//writeVTK2(t, currentfield, name.data(), w, h, startx, starty);
			}

#pragma omp barrier

#pragma omp master
			{
				Sleep(100);
				activeField = !activeField;
			}
		}
	}

	free(fields);
	free(otherfields);
}

int main(int c, char** v)
{
	srand(time(nullptr));

	int w = 0, h = 0;
	if (c > 1)
		w = atoi(v[1]); ///< read width
	if (c > 2)
		h = atoi(v[2]); ///< read height
	if (w <= 0)
		w = 30; ///< default width
	if (h <= 0)
		h = 30; ///< default height
	game(w, h);
}
