# Flex-MPI iterative matrix multiplication

This C programs perform distributed parallel matrix multiplication. Flex-MPI version (mat_x_adm) using MPI with support for dynamic process management via Flex-MPI (ADM framework), enabling processes to join or leave during execution while maintaining computation correctness and synchronization.

# Building ADM environment (Based on demo_outline.pdf)
1. Unzip PDP2023.zip
2. Execute command to extract:
   ```sh
   tar zxvf slurm_cluster_pdp_2023-03-01-4.tgz
4. Launch cluster:
   ```sh
   cd slurm_cluster_pdp
   ./launch-slurm-cluster.sh –n 3 –c 2
5. Copy Matrix multiplication programs
    ```sh
    docker cp MatrixMultiplication_Flex-MPI docker_image_id:/home/admin/tutorial_examples
6. Open shell in a cluster:
     ```sh
     docker exec -it slurm-node-1-1 /bin/bash
7. Change the Makefile and make:
   ```sh
   cd tutorial_examples
   nano Makefile
   //Add MatrixMultiplication_Flex-MPI to the folders list
   make

# Execute MPI version
1. Execute:
   ```sh
   cd /home/admin
   sbatch /home/admin/tutorial_examples/MatrixMultiplication_Flex-MPI/mat_x_mpi.sbatch

# Execute Flex-MPI version
1. Run slurm view script:
   ```sh
   icc_server
2. Execute:
   ```sh
   cd /home/admin
   sbatch /home/admin/tutorial_examples/MatrixMultiplication_Flex-MPI/mat_x_adm.sbatch

# Monitoring the execution out file
Opening another shell in any cluster this command shows the execution trace of the programs
  ```sh
  shared/view_results.sh
