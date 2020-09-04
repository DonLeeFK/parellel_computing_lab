#include<stdio.h>
#include<mpi.h>

int main(int argc, char* argv[]) {
    int Comm_rank, Comm_size;
    int i, n = 0;
    double sum, width, local, localpi, pi;
    double startTime = 0.0, endTime = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &Comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Comm_rank);

    if (Comm_rank == 0) {
        scanf("%d", &n);
        startTime = MPI_Wtime();
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    sum = 0.0;
    width = 1.0 / n;
    for (i = Comm_rank; i < n; i += Comm_size) {
        local = width * ((double)i + 0.5);
        sum += 4.0 / (1.0 + local * local);
    }
    localpi = width * sum;
    MPI_Reduce(&localpi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (Comm_rank == 0) {
        printf("pi = %.20f\n", pi);
        endTime = MPI_Wtime();
        printf("time = %f", endTime - startTime);
        //printf("rate = %f", (endTime - startTime)/0.000937);
    }

    MPI_Finalize();
    return 0;
}