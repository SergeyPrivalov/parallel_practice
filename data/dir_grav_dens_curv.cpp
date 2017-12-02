#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define GRAVY_CONST 6.67408
#define infilename "density2_100x100.dat"				
#define h1filename "h1_100x100.dat"		
#define h2filename "h2_100x100.dat"		
#define outfilename "field2_100x100.dat"

//#define outnoise "ds_noise_128x128.dat"				


bool ReadXYZFromFile (char * filename, double * x, double * y, double * z, int sizex, int sizey);
bool WriteXYZToFile(char * filename, double * x, double * y, double * z, int sizex, int sizey);
void addNoise(double *z, int size, double noiselevel);
void addNoise2(double *z, int size, double noiselevel, int noisenum);

int main() {
	srand (time(NULL));
	//double noiselevel=0.1;				

	int max_x=100;					
	int max_y = 100;					

	double dx,dy,xn,yn,x0,y0;
	double temp;

	double* xz,*yz,*z,*h1,*h2,*result;

	xz = (double *)malloc((max_x)*sizeof(double));
	yz = (double *)malloc((max_y)*sizeof(double));
	z = (double *)malloc((max_x*max_y)*sizeof(double));
	h1 = (double *)malloc((max_x*max_y) * sizeof(double));
	h2 = (double *)malloc((max_x*max_y) * sizeof(double));
	result = (double *)malloc((max_x*max_y)*sizeof(double));

	ReadXYZFromFile(infilename, xz,yz,z,max_x,max_y);
	ReadXYZFromFile(h1filename, xz, yz, h1, max_x, max_y);
	ReadXYZFromFile(h2filename, xz, yz, h2, max_x, max_y);

	x0 = xz[0];
	xn = xz[max_x-1];
	y0 = yz[0];
	yn = yz[max_y-1];
	dx = abs(xn-x0)/(max_x-1);
	dy = abs(yn-y0)/(max_y-1);


	double c = dx*dy*GRAVY_CONST;



	int zind;
	int resind = 0;
		for (int i = 0; i < max_y; i++) {
			printf("%i/%i\r",i+1,max_y);
			for (int j = 0; j < max_x; j++)
			{
			temp = 0.0;
				for (int ii = 0; ii < max_y; ii++) {
					for (int jj = 0; jj < max_x; jj++) {
						
				temp+= 
					(1.0/sqrt(
						(xz[jj]-xz[j])*(xz[jj]-xz[j])
						+(yz[ii]-yz[i])*(yz[ii]-yz[i])
						+h1[ii*max_x + jj] * h1[ii*max_x + jj]
					)
				-
				1.0/sqrt(
						(xz[jj]-xz[j])*(xz[jj]-xz[j])
						+(yz[ii]-yz[i])*(yz[ii]-yz[i])
						+h2[ii*max_x + jj] * h2[ii*max_x + jj]
					)
					)*z[ii*max_x+jj];
					;
				
	}}
				result[resind] = temp*c;
	resind++;
	}}

		printf("\n");
		WriteXYZToFile(outfilename, xz,yz,result,max_x,max_y);
	//	addNoise(result,int(max_x*max_y),0.15);
	//	WriteXYZToFile(outnoise, xz,yz,result,max_x,max_y);

	return 0;
}


bool ReadXYZFromFile (char * filename, double * x, double * y, double * z, int sizex, int sizey)
{
	FILE * fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Can not open file %s\n", filename);
		return false;
	}
	int k = 0;
	//printf ("%i, %i\n",sizex,sizey);
	for (int i = 0; i < sizey; i++)
		for (int j = 0; j < sizex; j++)
		{
			fscanf(fp, "%lf %lf %lf", x+j, y+i, z+k/*index(i,j,sizex)*/);
			k++;
		}
	fclose(fp);
	return true;
}

bool WriteXYZToFile (char * filename, double * x, double * y, double * z, int sizex, int sizey)
{
	FILE * fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("Can not open file %s\n", filename);
		return false;
	}
	int k = 0;
	//printf ("%i, %i\n",sizex,sizey);
	for (int i = 0; i < sizey; i++)
		for (int j = 0; j < sizex; j++)
		{
			fprintf(fp, "%5.20f %5.20f %5.20f\n", x[j], y[i], z[k]/*index(i,j,sizex)*/);
			k++;
		}
	fclose(fp);
	return true;
}

void addNoise(double *z, int size, double noiselevel)
{
	
	double min_n = 5555555.0;
	double max_n = -5555555.0;
	for (int i = 0; i < size; i++)
	{
		min_n > z[i] ? min_n = z[i] : 1;
		max_n < z[i] ? max_n = z[i] : 1;
	}
	printf("%f\n",max_n-min_n);
	double noiseval = (max_n-min_n)*noiselevel;
	 srand ( time(NULL) );

	for (int i = 0; i < size; i++)
	{
		if ( ((double)rand() / (double)RAND_MAX) > noiselevel) {
		z[i] += ((double)rand() / (double)RAND_MAX)*noiseval - noiseval/2;
		}
	}
}

void addNoise2(double *z, int size, double noiselevel, int noisenumber)
{
	for (int i = 0; i < size; i++) {
		z[i]=0.;
	}
	for (int i = 0; i < noisenumber; i++) {
		int pnt = ((double)rand() / (double)RAND_MAX)*size;
		z[pnt] = ((double)rand() / (double)RAND_MAX)*noiselevel;
	}
}




