#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define FILAS_A 5
#define COLUMNAS_A 3
#define FILAS_B 3
#define COLUMNAS_B 4
#define FILAS_C FILAS_A
#define COLUMNAS_C COLUMNAS_B

void generarMatrizAleatoria(int filas, int columnas, int matriz[filas][columnas]) {
    int i, j;
    
    for (i = 0; i < filas; i++) {
        for (j = 0; j < columnas; j++) {
            matriz[columnas*i+j] = rand() % 100;  // Generar números aleatorios entre 0 y 99
        }
    }
}
/*
void multiplicarMatrices(int filas_A, int columnas_A, int columnas_B, int matriz_A[filas_A][columnas_A], int matriz_B[columnas_A][columnas_B], int matriz_C[filas_A][columnas_B]) {
    int i, j, k;

    for (i = 0; i < filas_A; i++) {
        for (j = 0; j < columnas_B; j++) {
            matriz_C[i][j] = 0;
            for (k = 0; k < columnas_A; k++) {
                matriz_C[i][j] += matriz_A[i][k] * matriz_B[k][j];
            }
        }
    }
}

void multiplicarMatrices2(int filas_A, int columnas_A, int columnas_B, int matriz_A[filas_A][columnas_A], int matriz_B[columnas_A][columnas_B], int matriz_C[filas_A][columnas_B]) {
    int i, j, k, a, b;

    for (int i = inicio; i < fin; i++)
    {
        matriz_C[i] = 0;
        a = i/filas_A;
        b = i%columnas_B;
        for (int j = 0; j < columnas_A; j++)
        {
            matriz_C[i] += matriz_A[columnas_A*a+j]*matriz_B[columnas_B*j+b];
        }
    }
}
*/
void imprimirMatriz(int filas, int columnas, int matriz[filas][columnas]) {
    int i, j;

    for (i = 0; i < filas; i++) {
        for (j = 0; j < columnas; j++) {
            printf("%3d ", matriz[columnas*i+j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int matriz_A[FILAS_A][COLUMNAS_A];
    int matriz_B[FILAS_B][COLUMNAS_B];
    int matriz_C[FILAS_C][COLUMNAS_C];

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
/*
    if (size != 1 && size != FILAS_C) {
        if (rank == 0) {
            printf("El número de procesos debe ser 1 o igual al número de filas de la matriz resultante.\n");
        }
        MPI_Finalize();
        return 0;
    }
*/
    if (rank == 0) {
        srand(time(NULL));  // Inicializar la semilla aleatoria
        generarMatrizAleatoria(FILAS_A, COLUMNAS_A, matriz_A);
        generarMatrizAleatoria(FILAS_B, COLUMNAS_B, matriz_B);
    }

    MPI_Bcast(matriz_A, FILAS_A * COLUMNAS_A, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(matriz_B, FILAS_B * COLUMNAS_B, MPI_INT, 0, MPI_COMM_WORLD);

    int inicio = ((FILAS_C*COLUMNAS_C) / size) * rank;
    int fin = inicio+(FILAS_C*COLUMNAS_C)/size;
    int a, b;

    for (int i = inicio; i < fin; i++)
    {
        matriz_C[i] = 0;
        a = i/FILAS_A;
        b = i%COLUMNAS_B;
        for (int j = 0; j < COLUMNAS_A; j++)
        {
            matriz_C[i] += matriz_A[COLUMNAS_A*a+j]*matriz_B[COLUMNAS_B*j+b];
        }
    }

/*
    int submatriz_C[fin - inicio][COLUMNAS_C];

    multiplicarMatrices(fin - inicio, COLUMNAS_A, COLUMNAS_B, matriz_A + inicio, matriz_B, submatriz_C);
*/
    int *buffer = NULL;
    if (rank == 0) {
        buffer = (int *)malloc(FILAS_C * COLUMNAS_C * sizeof(int));
    }

    MPI_Gather(resultado + inicio, fin - inicio, MPI_INT, buffer, fin - inicio, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Matriz A:\n");
        imprimirMatriz(FILAS_A, COLUMNAS_A, matriz_A);

        printf("\nMatriz B:\n");
        imprimirMatriz(FILAS_B, COLUMNAS_B, matriz_B);

        printf("\nResultado de la multiplicación:\n");
        imprimirMatriz(FILAS_C, COLUMNAS_C, buffer);

        free(buffer);
    }

    MPI_Finalize();
    return 0;
}
