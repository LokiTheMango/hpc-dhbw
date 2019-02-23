#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <string>
#include <fstream>

#define min(x, y) (x > y ? y : x)

#define calcIndex(width, x, y) ((y) * (width) + (x))

long TimeSteps = 100;

void writeVTK2(long timestep, double *data, const char *prefix, int w, int h, int offsetX, int offsetY)
{
  char filename[2048];
  int x, y;

  float deltax = 1.0;
  long nxy = w * h * sizeof(float);

  snprintf(filename, sizeof(filename), "%s-%05ld%s", prefix, timestep, ".vti");
  FILE *fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
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
      float value = data[calcIndex(h, x, y)];
      fwrite((unsigned char *)&value, sizeof(float), 1, fp);
    }
  }

  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void show(double *currentfield, int w, int h)
{
  printf("\033[H");
  int x, y;
  for (y = 0; y < h; y++)
  {
    for (x = 0; x < w; x++)
      printf(currentfield[calcIndex(w, x, y)] ? "\033[07m  \033[m" : "  ");
    //printf("\033[E");
    printf("\n");
  }
  fflush(stdout);
}

void evolve(double *currentfield, double *newfield, int w, int h, int start_x, int start_y, int end_x, int end_y)
{
  int x, y;
  for (y = start_y; y < min(end_y, h); y++)
  {
    for (x = start_x; x < min(end_x, w); x++)
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

void filling(double *currentfield, int w, int h)
{
  int i;
  for (i = 0; i < h * w; i++)
  {
    currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
  } return;

  std::ifstream im("im.txt");
  im.seekg(0, std::ios::end);
  auto len = im.tellg();
  im.seekg(0, std::ios::beg);

  
  std::string field;
  field.resize(len);
  im.read((char*)field.data(), len);
  im.close();

  int pos = 0;
  for(int i = 0; i < h * w; ++i) {
    char x = field[pos++ % field.size()];
    if(x == '0') {
      currentfield[i] = 0.0;
    }
   else  if(x == '1') {
      currentfield[i] = 1.0;
    }
    else {
      --i;
    }
  }
}

void game(int w, int h)
{
  double *currentfield = (double *)calloc(w * h, sizeof(double));
  double *newfield = (double *)calloc(w * h, sizeof(double));

  //printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

  filling(currentfield, w, h);
  long t;
  for (t = 0; ; t++)
  {
    show(currentfield, w, h);


#pragma omp parallel
    {
      int threadnums = omp_get_num_threads();
      int threadnum = omp_get_thread_num();
      int startx=15*(threadnum<2);
      int starty=15*(threadnum%2);
      int endx=startx+15;
      int endy=starty+15;
      evolve(currentfield, newfield, w, h, startx, starty, endx, endy);

      printf("%ld timestep\n", t);

      std::string name("gol-" + std::to_string(omp_get_thread_num()));
      //writeVTK2(t, currentfield, name.data(), w, h, startx, starty);
    }
    usleep(100000);

    //SWAP
    double *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v)
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
