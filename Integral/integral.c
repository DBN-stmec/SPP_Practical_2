//  Author:
//
//    Yannic Fischler
//    Modfied for WS21 by Sebastian Kreutzer

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define UPPER_LIMIT 1
#define LOWER_LIMIT -1
#define INTEGRAL(x) (2/(1+x*x))

int world_size;
int world_rank;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // Store number of processes in world_size
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Store current rank in world_rank
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc != 2) {
        printf("Add number of sampling points as parameter\n");
        return 1;
    }
    int numberOfPoints = atoi(argv[1]);
    // Make sure that numberOfPoints is valid
    while (numberOfPoints < world_size) {
        if(world_rank == 0){
            printf("Too few Datapoints.\n Please enter a higher number than: %d \n",world_size);
            scanf("%d",&numberOfPoints);
        
        for(int i = 1; i<world_size;i++){
            MPI_Send(&numberOfPoints,1,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        }
        if(world_rank != 0){
            MPI_Recv(&numberOfPoints,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        //MPI_Bcast(&numberOfPoints , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
        if (numberOfPoints >= world_size){
            break;
        }
    }
    if (world_rank == 0) {
        printf("Running with %d processes\n", world_size);
    }
    double result = 0;
    double totalSum = 0;
    double sum;
    // Calculate the number of samples for each process
    int localSamples = numberOfPoints / world_size;

    for (int i = world_rank * localSamples; i <= world_rank * localSamples + localSamples -1; i++) {
        // Generate a random point in the interval [-1, 1]
        double x = LOWER_LIMIT + (UPPER_LIMIT-LOWER_LIMIT) * (double)rand() / (double)RAND_MAX;
        sum += INTEGRAL(x); // Evaluate the function at the point
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank != world_size-1 && world_rank != 0){
        MPI_Recv(&totalSum, 1, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        totalSum += sum;
        MPI_Send(&totalSum, 1, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
    }
    // Send the local sum to the first process
    if (world_rank == world_size-1) { // Slave thread works
        MPI_Send(&sum, 1, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);              
    } 
     // Master thread works
        // The first process receives the local sums from all other processes
        //result = totalSum;
        //for (int i = 1; i < world_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        //double sum;

        //}
        //MPI_Barrier(MPI_COMM_WORLD);
        // Calculate and print the value of the integral
        if(world_rank ==0){
            MPI_Recv(&totalSum, 1, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalSum += sum;
            result = (totalSum / numberOfPoints) * (UPPER_LIMIT-LOWER_LIMIT);  // The interval has length UPPER_LIMIT-LOWER_LIMIT
            printf("%f\n", result);  
        }

    

    MPI_Finalize();
    return 0;

}
