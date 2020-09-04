#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define num_threads 8
#define N 100000

int main(int argc, char* argv[])
{
    int i;
    double local, tmp, pi, w;
    double startTime, end;
    omp_set_num_threads(num_threads);
    w = 1.0 / N;
    pi = 0;
    startTime = omp_get_wtime();
    #pragma omp parallel for default(shared) private(i, local, tmp) reduction(+ : pi)
    for (i = 0; i < N; i++)
    {
        local = (i + 0.5) * w;
        tmp = 4.0 / (1.0 + local * local);
        pi += tmp;
    }
    end = omp_get_wtime();
    pi = pi * w;
    printf("The value of pi is %.12lf\n", pi);
    printf("the time cost is %lfs\n", end - startTime);
    printf("%lf\n", end - startTime);
    return 0;
}