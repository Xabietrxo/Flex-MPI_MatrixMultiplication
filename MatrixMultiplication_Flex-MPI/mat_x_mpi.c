#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROWS1 17
#define COLS1 16
#define ROWS2 16
#define COLS2 15

void generateRandomMatrix(int *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[cols * i + j] = rand() % 100;
        }
    }
}

void printMatrix(int *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%3d ", matrix[cols * i + j]);
        }
        printf("\n");
    }
}

void mulMatrix(int *matrix1, int *matrix2, int *resultMatrix,
               int rows1, int cols1, int cols2) {
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            resultMatrix[i * cols2 + j] = 0;
            for (int k = 0; k < cols1; k++) {
                resultMatrix[i * cols2 + j] += matrix1[i * cols1 + k] * matrix2[k * cols2 + j];
            }
        }
    }
}

void compareMatrix(int *matrix1, int *matrix2, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%s |", (matrix1[cols * i + j] == matrix2[cols * i + j]) ? "true" : "false");
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    int it=0;
    int itmax=ROWS1*COLS2;

    int *matrix1 = malloc(ROWS1 * COLS1 * sizeof(int));
    int *matrix2 = malloc(ROWS2 * COLS2 * sizeof(int));
    int *result  = malloc(ROWS1 * COLS2 * sizeof(int));

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        srand(time(NULL));  // Initiate random seed
        generateRandomMatrix(matrix1, ROWS1, COLS1);
        generateRandomMatrix(matrix2, ROWS2, COLS2);
    }

    MPI_Bcast(matrix1, ROWS1 * COLS1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(matrix2, ROWS2 * COLS2, MPI_INT, 0, MPI_COMM_WORLD);

    int local_count = 0;

    for (; it < itmax; it++)
    {
        if (it%size == rank)
        {
            result[it] = 0;
            for (int k = 0; k < COLS1; k++) {
                result[it] += matrix1[(it/COLS2) * COLS1 + k] * matrix2[k * COLS2 + (it%COLS2)];
            }
            local_count++;
        }
    }

    MPI_Datatype intercalate_type;
    MPI_Type_vector(local_count, 1, size, MPI_INT, &intercalate_type);
    MPI_Type_commit(&intercalate_type);

    if (rank!=0)
    {
        MPI_Send(&result[rank], 1, intercalate_type, 0, 0, MPI_COMM_WORLD);
    } else{
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(&result[i], 1, intercalate_type, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0) {

        int *resultMatrix = malloc(ROWS1 * COLS2 * sizeof(int));
        mulMatrix(matrix1, matrix2, resultMatrix, ROWS1, COLS1, COLS2);

        printf("Matrix 1:\n");
        printMatrix(matrix1, ROWS1, COLS1);

        printf("\nMatrix 2:\n");
        printMatrix(matrix2, ROWS2, COLS2);

        printf("\nMultiplication result:\n");
        printMatrix(result, ROWS1, COLS2);

        printf("\nCompare:\n");
        compareMatrix(result, resultMatrix, ROWS1, COLS2);

        free(result);
        free(resultMatrix);
    }

    MPI_Finalize();
    return 0;
}
