#include <iostream>
#include <omp.h>

void FillMatrix(int** matrix, int n, int m, int number = 1)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			matrix[i][j] = number;
		}
	}
}

bool IsResultCorrect(int** a, int** b, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (a[i][j] != b[i][j])
				return 0;
		}
	}
	return 1;
}

int main()
{
	long n = 2000;
	long n1 = n;
	long m1 = n;
	long n2 = n;
	long m2 = n;

	int** a = new int* [n1]; 
	for (int i = 0; i < n1; i++)
		a[i] = new int[m1];

	int** b = new int* [n2];
	for (int i = 0; i < n2; i++)
		b[i] = new int[m2];

	int** c = new int* [n1];
	for (int i = 0; i < n1; i++)
		c[i] = new int[m2];

	int** c_parallel = new int* [n1];
	for (int i = 0; i < n1; i++)
		c_parallel[i] = new int[m2];

	FillMatrix(a, n1, m1, 2);
	FillMatrix(b, n2, m2, 3);
	FillMatrix(c, n1, m2, 0);
	FillMatrix(c_parallel, n1, m2, 0);


	double time_start = omp_get_wtime();      // POSLEDOVATEL'NIY
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < m2; j++) {
			for (int k = 0; k < m1; k++) {
				c[i][j] += (a[i][k] * b[k][j]);
			}
		}
	}
	double time_end = omp_get_wtime();
	double seq_time_result = 1000 * (time_end - time_start);
	std::cout << "Sequential code time action is " << seq_time_result << " ms\n";



	omp_set_num_threads(2);   // 2 THREDA
	int i, j, k;
	time_start = omp_get_wtime();
#pragma omp parallel for shared(a, b, c_parallel) private(i, j, k)
	for (i = 0; i < n1; i++) {
		for (j = 0; j < m2; j++) {
			for (k = 0; k < m1; k++) {
				c_parallel[i][j] += (a[i][k] * b[k][j]);
			}
		}
	}
	time_end = omp_get_wtime();
	double par_time_result = 1000 * (time_end - time_start);
	std::cout << "Parallel 2 threads code time action is " << par_time_result << " ms\n";
	if (!IsResultCorrect(c, c_parallel, n1, m2))
		std::cout << "Matrix is not correct\n";
	


	FillMatrix(c_parallel, n1, m2, 0);      
	omp_set_num_threads(4);                  // 4 THREDA
	time_start = omp_get_wtime();
#pragma omp parallel for shared(a, b, c_parallel) private(i, j, k)
	for (i = 0; i < n1; i++) {
		for (j = 0; j < m2; j++) {
			for (k = 0; k < m1; k++) {
				c_parallel[i][j] += (a[i][k] * b[k][j]);
			}
		}
	}
	time_end = omp_get_wtime();
	par_time_result = 1000 * (time_end - time_start);
	std::cout << "Parallel 4 threads code time action is " << par_time_result << " ms\n";
	if (!IsResultCorrect(c, c_parallel, n1, m2))
		std::cout << "Matrix is not correct\n";


	for (int i = 0; i < n1; i++)
		delete[]a[i];
	delete[]a;
	for (int i = 0; i < n2; i++)
		delete[]b[i];
	delete[]b;
	for (int i = 0; i < n1; i++)
		delete[]c[i];
	delete[]c;
	for (int i = 0; i < n1; i++)
		delete[]c_parallel[i];
	delete[]c_parallel;
	system("pause");
	return 0;
}