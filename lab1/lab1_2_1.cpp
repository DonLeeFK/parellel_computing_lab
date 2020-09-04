#include<stdio.h>
#include<omp.h>
#include<math.h>

#define num_threads 8
#define N 1000

void main()
{
    int m, k, i, count = 0;
    double startTime, endTime, time;
    omp_set_num_threads(num_threads);
    startTime = omp_get_wtime();
#pragma omp parallel for default(shared) reduction(+:count) private(i,k,m)
    for (m = 2; m <= N; m = m + 1)
    {
        k = sqrt(m);
        for (i = 2; i <= k; i++)
        {
            if (m % i == 0) break;
        }
        if (i > k)
        {
#pragma omp atomic              
            count++;
        }
    }
    endTime = omp_get_wtime();
    time = (endTime - startTime);
    printf("count = %d\ntime = %f\n", count, time);
}