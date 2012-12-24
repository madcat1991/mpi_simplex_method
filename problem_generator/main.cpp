#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

bool parseArgs(int argc, char *argv[], int &n, int &m, char *filename);
void printHelp();

int main(int argc, char* argv[])
{
	const int GEN_RANGE_MAX = 100;
	const int GEN_RANGE_MIN = -100;

	int n = 0;     // Number of variables
	int m = 0;     // Number of constraints
	int count = 0; // Number of constraint's coefficients  
	char filename[FILENAME_MAX];
	FILE *f;
	double *random; // Constraint's coefficients
	double sum_sqrs = 0;
	double scale = 0;
	int i, j;
	double coef;
	double b;

	if (!parseArgs(argc, argv, n, m, filename))
	{
		return 1;
	}

	f = fopen(filename, "w");

	if (f == NULL)
	{
		printf("Cannot create file: %s", filename);
		return 2;
	}

	count = n * m;

	random = (double *) malloc(sizeof(double) * count);
	
	/*srand((unsigned int) time(NULL));*/
	srand(234);

	for (i = 0; i < count; i++)
	{
		random[i] = (double)rand() / (GEN_RANGE_MAX + 1) * (GEN_RANGE_MAX - GEN_RANGE_MIN) + GEN_RANGE_MIN;
		sum_sqrs += random[i] * random[i];
	}

	scale = sqrt(sum_sqrs) / 100;

	for (i = 0; i < count; i++)
	{
    	    random[i] /= scale;
	}

	fprintf(f, "%d\n", n);
	fprintf(f, "%d\n", m);

    // Coefficients of variables in function
	for (i = 0; i < n; i++)
	{
		coef = (double)rand() / (GEN_RANGE_MAX + 1) * (GEN_RANGE_MAX - GEN_RANGE_MIN) + GEN_RANGE_MIN;
    		fprintf(f, "%lf\n", coef);                   
	}

	// Hyperplane: A1*(x1-A1) + A2*(x2-A2) + ... + An*(xn - An) = 0
	// then A1*x1 + A2*x2 + ... + An*xn = sum(Ai*Ai), where 0<i<=n
    
	// Coefficients of constraints
	for (i = 0; i < m; i++)
	{
    	    b = 0;
        // Coefficients of variables in constraint
    	    for (j = 0; j < n; j++)
    	    {
		fprintf(f, "%lf\n", random[i + j]);
        	b += random[i + j] * random[i + j];
    	    }
		
		// Free coefficient of constraint
		fprintf(f, "%lf\n", b);
	}
}


bool parseArgs(int argc, char *argv[], int &n, int &m, char *filename)
{
	const char *DEFAULT_FILENAME = "input.txt";
	if (argc < 2)
	{
		printHelp();
		return false;
	}

	n = atoi(argv[1]);
	
	if (argc > 2)
		m = atoi(argv[2]);
	else
		m = n;

	if (argc > 3)
		strcpy(filename, argv[3]);
	else
		strcpy(filename, DEFAULT_FILENAME);

	if (n <= 0 || m <= 0)
	{
		printf("Incorrect number!\n");
		printHelp();
		return false;
	}

	return true;
}

void printHelp()
{
	printf("Usage: generate n [m] [filename]\n");
	printf("Parameters: \n");
	printf("\tn\t\tNumber of variables.\n");
	printf("\tm\t\tNumber of constraints.\n");
	printf("\tfilename\tOutput filename.\n");
	printf("\n");
	printf("If m is not specified then m = n.\n");
	printf("If filename is not specified then output filename is \"input.txt\".\n");
	printf("\n");
	printf("Structure of generated file: \n");
	printf("n\t\t  - number of variables\n");
	printf("m\t\t  - number of constraints\n");
	printf("c[0]\t\t  - coefficients of variables in function\n");
	printf("..\t\t  - \n");
	printf("c[n-1]\t\t  - \n");
	printf("A[0, 0]\t\t  - coefficients of variables in constraint #1\n");
	printf("..\t\t  - \n");
	printf("A[m-1, n-1]\t  - \n");
	printf("b[0]\t\t  - free coefficient in constraint #1\n");
	printf("..\n");
	printf("A[m-1][0]\t  - coefficients of variables in constraint #m\n");
	printf("..\t\t  - \n");
	printf("A[m-1, n-1]\t  - \n");
	printf("b[m-1]\t\t  - free coefficient in constraint #m\n");
	printf("\n");
	printf("\n");
}
