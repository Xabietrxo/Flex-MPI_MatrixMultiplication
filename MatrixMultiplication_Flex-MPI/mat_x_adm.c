#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <empi.h>
#include <stdbool.h>


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
            printf("%s ", (matrix1[cols * i + j] == matrix2[cols * i + j]) ? "true" : "false");
        }
        printf("\n");
    }
}


void evaluate(int it, int itmax, bool *sync, bool *send, int *procs_hint, int *excl_nodes_hint);
void isWorldSizeChanged(int world_size, int *last_world_size, int *saved_last_world_size, int *matrix1, int *matrix2, int dim1, int dim2, bool sync, bool *send);
void compute(int it, int last_it, int world_size, int world_rank, int *result, int *matrix1, int *matrix2, int cols1, int cols2, int *ec);
void gather(bool *send, bool *sync, int last_world_size, int saved_last_world_size, int world_size, int world_rank, int *ec, int *result, int it, int *last_it);

int main(int argc, char *argv[]) {

    int it=0;
    int itmax=ROWS1*COLS2;

    int *matrix1 = malloc(ROWS1 * COLS1 * sizeof(int));
    int *matrix2 = malloc(ROWS2 * COLS2 * sizeof(int));
    int *result  = malloc(ROWS1 * COLS2 * sizeof(int));

    int world_rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(ADM_COMM_WORLD, &world_rank);
    MPI_Comm_size(ADM_COMM_WORLD, &world_size);

    int dim1 = ROWS1*COLS1;
    int dim2 = ROWS2*COLS2;

    if (world_rank == 0) {
        srand(time(NULL));
        generateRandomMatrix(matrix1, ROWS1, COLS1);
        generateRandomMatrix(matrix2, ROWS2, COLS2);
    }

    MPI_Bcast(matrix1, dim1, MPI_INT, 0, ADM_COMM_WORLD);    //Broadcast both matrix
    MPI_Bcast(matrix2, dim2, MPI_INT, 0, ADM_COMM_WORLD);

    // get process type
    int proctype;
    ADM_GetSysAttributesInt ("ADM_GLOBAL_PROCESS_TYPE", &proctype);

    // if process is native
    if (proctype == ADM_NATIVE) {
        
        printf ("Rank(%d/%d): Process native\n", world_rank, world_size);
        
    // if process is spawned
    } else {
        printf ("Rank(%d/%d): Process spawned\n", world_rank, world_size);
    }
    
    /* set max number of iterations */
    ADM_RegisterSysAttributesInt ("ADM_GLOBAL_MAX_ITERATION", &itmax);

    /* get actual iteration for new added processes*/
    ADM_GetSysAttributesInt ("ADM_GLOBAL_ITERATION", &it);

    // init last world size, flags, counters...
    int last_world_size = world_size;
    int saved_last_world_size = world_size;
    int last_it = it;
    bool send = false;
    bool sync = false;
    int ec = 0;

    // start loop
    for (; it < itmax; it++)
    {
        int procs_hint = 0;
        int excl_nodes_hint = 0;

        evaluate(it, itmax, &sync, &send, &procs_hint, &excl_nodes_hint);

        isWorldSizeChanged(world_size, &last_world_size, &saved_last_world_size, matrix1, matrix2, dim1, dim2, sync, &send);

        gather(&send, &sync, last_world_size, saved_last_world_size, world_size, world_rank, &ec, result, it, &last_it);
        
        
        /* start malelability region */
        ADM_MalleableRegion (ADM_SERVICE_START);

        compute(it, last_it, world_size, world_rank, result, matrix1, matrix2, COLS1, COLS2, &ec);
        
        // update the iteration value
        ADM_RegisterSysAttributesInt ("ADM_GLOBAL_ITERATION", &it);        
        
        // ending malleable region
        int status;
        status = ADM_MalleableRegion (ADM_SERVICE_STOP);

        
        if (status == ADM_ACTIVE) {
            // updata rank and size
            MPI_Comm_rank(ADM_COMM_WORLD, &world_rank);
            MPI_Comm_size(ADM_COMM_WORLD, &world_size);
        } else {
            printf("Rank (%d/%d): Iteration:= %d, Break\n", world_rank, world_size, it);
            // end the process
            break;
        }
    }

    printf("Rank (%d/%d): End of loop \n", world_rank, world_size);

    MPI_Finalize();

    if (world_rank == 0) {
        int *resultMatrix = malloc(ROWS1 * COLS2 * sizeof(int));
        mulMatrix(matrix1, matrix2, resultMatrix, ROWS1, COLS1, COLS2);

        printf("Matrix 1:\n");
        printMatrix(matrix1, ROWS1, COLS1);

        printf("\nMatrix 2:\n");
        printMatrix(matrix2, ROWS2, COLS2);

        printf("\nMultiplication result (non-parallel):\n");
        printMatrix(resultMatrix, ROWS1, COLS2);

        printf("\nMultiplication result (parallel):\n");
        printMatrix(result, ROWS1, COLS2);

        printf("\nCompare:\n");
        compareMatrix(result, resultMatrix, ROWS1, COLS2);

        free(resultMatrix);
    }

    
    return 0;
}

void evaluate(int it, int itmax, bool *sync, bool *send, int *procs_hint, int *excl_nodes_hint){

    if ( (it == 2*(itmax/10)) || (it == 4*(itmax/10)) ){
        *procs_hint = 2;
        *excl_nodes_hint = 0;
        *sync = true;
        ADM_RegisterSysAttributesInt ("ADM_GLOBAL_HINT_NUM_PROCESS", procs_hint);
        ADM_RegisterSysAttributesInt ("ADM_GLOBAL_HINT_EXCL_NODES", excl_nodes_hint);
        printf("Rank (%d/%d): Iteration:= %d, procs_hint=%d, excl_nodes_hint=%d\n", world_rank, world_size, it, procs_hint, excl_nodes_hint);
    } else if ( (it == 6*(itmax/10)) || (it == 8*(itmax/10)) ){
        *procs_hint = -2;
        *excl_nodes_hint = 0;
        *send = true;
        ADM_RegisterSysAttributesInt ("ADM_GLOBAL_HINT_NUM_PROCESS", procs_hint);
        ADM_RegisterSysAttributesInt ("ADM_GLOBAL_HINT_EXCL_NODES", excl_nodes_hint);
        printf("Rank (%d/%d): Iteration:= %d, procs_hint=%d, excl_nodes_hint=%d\n", world_rank, world_size, it, procs_hint, excl_nodes_hint);
    }else if (it == (itmax-1)){
        *send = true;
    }
}

void isWorldSizeChanged(int world_size, int *last_world_size, int *saved_last_world_size, int *matrix1, int *matrix2, int dim1, int dim2, bool sync, bool *send){
    if (*last_world_size != world_size) {
        if (*last_world_size<world_size)
        {
            MPI_Bcast(matrix1, dim1, MPI_INT, 0, ADM_COMM_WORLD);    //Broadcast both matrix to spawned processes
            MPI_Bcast(matrix2, dim2, MPI_INT, 0, ADM_COMM_WORLD);
        }
        *saved_last_world_size = *last_world_size;
        *last_world_size = world_size;
        if (sync == true)
        {
            *send = true;
        }
    }
}

void compute(int it, int last_it, int world_size, int world_rank, int *result, int *matrix1, int *matrix2, int cols1, int cols2, int *ec){
    if ((it-last_it)%world_size == world_rank)        //Divide workload in processes
    {
        result[it] = 0;
        for (int k = 0; k < cols1; k++) {
            result[it] += matrix1[(it/cols2) * cols1 + k] * matrix2[k * cols2 + (it%cols2)];
        }
        *ec = *ec + 1;
        printf("Rank (%d/%d): Iteration(last_it:%d):= %d, [Row: %d , Col: %d] = %d\n", world_rank, world_size, last_it, it, (it/cols2), (it%cols2), result[it]);
        sleep(1);
    }
}

void gather(bool *send, bool *sync, int last_world_size, int saved_last_world_size, int world_size, int world_rank, int *ec, int *result, int it, int *last_it){
    MPI_Datatype intercalate_type;

    if (*send == true)
    {
        if (*sync == false)
        {
            MPI_Type_vector(*ec, 1, last_world_size, MPI_INT, &intercalate_type);
            MPI_Type_commit(&intercalate_type);
            if (world_rank!=0)
            {
                MPI_Send(&result[*last_it+world_rank], 1, intercalate_type, 0, 0, ADM_COMM_WORLD);
                *ec = 0;
                printf("Rank (%d/%d): Iteration:= %d, Data sent\n", world_rank, world_size, it);
            }
            if (world_rank==0)
            {
                for (int i = 1; i < last_world_size; i++)
                {
                    MPI_Recv(&result[*last_it+i], 1, intercalate_type, i, 0, ADM_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }else{
            MPI_Type_vector(*ec, 1, saved_last_world_size, MPI_INT, &intercalate_type);
            MPI_Type_commit(&intercalate_type);
            if (world_rank!=0)
            {
                MPI_Send(&result[*last_it+world_rank], 1, intercalate_type, 0, 0, ADM_COMM_WORLD);
                *ec = 0;
                printf("Rank (%d/%d): Iteration:= %d, Data sent\n", world_rank, world_size, it);
            }
            if (world_rank==0)
            {
                for (int i = 1; i < saved_last_world_size; i++)
                {
                    MPI_Recv(&result[*last_it+i], 1, intercalate_type, i, 0, ADM_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            *sync = false;
        }
        *last_it = it;
        *send = false;
    }

    //MPI_Type_free(&intercalate_type);
}