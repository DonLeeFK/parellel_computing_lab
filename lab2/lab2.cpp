#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define v_max 20
#define car 1000000
#define cycle 300
#define p 0.5

int main(int argc, char* argv[]) {
    int size, rank, * v;
    long* position, distance;
    double randNum;
    double starttime, endtime;
    position = (long*)malloc(car * sizeof(long));
    v = (int*)malloc(car * sizeof(int));

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    starttime = MPI_Wtime();
    for (int j = rank * car / size; j < rank * car / size + car / size; j++) {
        position[j] = j;
        v[j] = 0;
    }
    if (size > 1) {
        for (int j = 0; j < size; j++) {
            MPI_Bcast(&position[j * car / size], car / size, MPI_LONG, j, MPI_COMM_WORLD);
            MPI_Bcast(&v[j * car / size], car / size, MPI_INT, j, MPI_COMM_WORLD);
        }
    }
    for (int i = 0; i < cycle; i++) {
        for (int j = rank * car / size; j < rank * car / size + car / size; j++) {
            if (j == car - 1 && v[j] < v_max) {
                v[j] ++;
            }
            else if (j != car - 1) {
                if (v[j] < v_max) {
                    v[j] ++;
                }
                distance = position[j + 1] - position[j];
                if (distance <= v[j]) {
                    v[j] = distance - 1;
                }
            }
            randNum = rand() / (double)RAND_MAX;
            if (randNum < p) {
                if (v[j] > 0) {
                    v[j]--;
                }
            }
            if (j == car - 1)
            {
                position[j] += v[j];
            }
            else
            {
                if (position[j] + v[j] < position[j + 1])
                {
                    position[j] += v[j];
                }
                else
                    position[j] = position[j + 1] - 1;
            }
        }
        if (size > 1) {
            for (int j = 0; j < size; j++) {
                MPI_Bcast(&position[j * car / size], car / size, MPI_LONG, j, MPI_COMM_WORLD);
                MPI_Bcast(&v[j * car / size], car / size, MPI_INT, j, MPI_COMM_WORLD);
            }
        }
    }
    endtime = MPI_Wtime();
    if (rank == 0) {
        for (int i = 0; i < car; i++) {
            printf("car:%d position:%ld v:%d\n", i, position[i], v[i]);
        }
        printf("time spent: %fs\n", endtime-starttime);
    }
    MPI_Finalize();
    return 0;
}