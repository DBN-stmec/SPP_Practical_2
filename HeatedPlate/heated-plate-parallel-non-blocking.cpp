#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>

#include <mpi.h>

using namespace std;

// Simple helper class to measure time.
struct Timer
{
  using clock = std::chrono::steady_clock;

  clock::time_point startTime;

  Timer() : startTime(clock::now()) {}

  double getSeconds()
  {
    return std::chrono::duration_cast<std::chrono::duration<float>>(clock::now() - startTime).count();
  }
};

// Verify result
bool verify(double *result, double *ref, long size)
{
  for (int i = 0; i < size; i++)
  {
    
    if (fabs(result[i] - ref[i]) >= 1e-3)      
      return false;
  }
  return true;
}

// Write reference file
bool writeToFile(std::string file, double *vals, long size)
{
  std::fstream fout(file, std::ios_base::out);
  if (!fout.is_open())
  {
    std::cout << "Unable to open file" << std::endl;
    return false;
  }
  for (int i = 0; i < size; i++)
  {
    fout << vals[i] << " ";
  }
  return true;
}

// Read reference file
bool readFromFile(std::string file, double *vals, long size)
{
  std::fstream fin(file, std::ios_base::in);
  long idx = 0;
  double val;
  while (fin >> val)
  {
    if (idx >= size)
    {
      std::cout << "Too many values in reference file" << std::endl;
      return false;
    }
    vals[idx++] = val;
  }
  if (size != idx)
  {
    std::cout << "Too few values in reference file" << std::endl;
    return false;
  }
  return true;
}

int world_size;
int world_rank;

int main(int argc, char *argv[])

//
//  Purpose:
//
//    MAIN is the main program for HEATED_PLATE_MPI.
//
//  Discussion:
//
//    This code solves the steady state heat equation on a rectangular region.
//
//    The sequential version of this program needs approximately
//    18/epsilon iterations to complete.
//
//
//    The physical region, and the boundary conditions, are suggested
//    by this diagram;
//
//                   W = 0
//             +------------------+
//             |                  |
//    W = 100  |                  | W = 100
//             |                  |
//             +------------------+
//                   W = 100
//
//    The region is covered with a grid of M by N nodes, and an M by N
//    array W is used to record the temperature.  The correspondence between
//    array indices and locations in the region is suggested by giving the
//    indices of the four corners:
//
//                  I = 0
//          [0][0]-------------[0][N-1]
//             |                  |
//      J = 0  |                  |  J = N-1
//             |                  |
//        [M-1][0]-----------[M-1][N-1]
//                  I = M-1
//
//    The steady state solution to the discrete heat equation satisfies the
//    following condition at an interior grid point:
//
//      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    where "Central" is the index of the grid point, "North" is the index
//    of its immediate neighbor to the "north", and so on.
//
//    Given an approximate solution of the steady state heat equation, a
//    "better" solution is given by replacing each interior point by the
//    average of its 4 neighbors - in other words, by using the condition
//    as an ASSIGNMENT statement:
//
//      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    If this process is repeated often enough, the difference between
//    successive estimates of the solution will go to zero.
//
//    This program carries out such an iteration, using a tolerance specified by
//    the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2011
//
//  Author:
//
//    Original C version by Michael Quinn.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
//  Local parameters:
//
//    Local, double DIFF, the norm of the change in the solution from one
//    iteration to the next.
//
//    Local, double MEAN, the average of the boundary values, used to initialize
//    the values of the solution in the interior.
//
//    Local, double U[M][N], the solution at the previous iteration.
//
//    Local, double W[M][N], the solution computed at the latest iteration.
//
{
#define M 500
#define N 500

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Request request;
  double diff;
  double partialDiff;
  double epsilon = 0.001;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;

  // TODO: Each process may only allocate a part of the matrix (DONE)
  // you may want to allocate the matrix wit malloc depening on the number of
  // lines each process may get
  int modulo = M%world_size;          //if number of rows cannot be divided "cleanly" by number of processes remainder(modulo) is added to number of rows in master process
  int subRows = (M-modulo)/world_size;
  if(world_rank==0){
    subRows += modulo;
  }
   // Number of rows(M) per process (world_rank)
  //double (*u)[N];
  double u[subRows][N];
  double w[subRows][N];
  //double** u = (double **)malloc(sizeof(double*) *subRows);
  //double **u =(double **)malloc(sizeof(double *) * N * subRows);
  //double** w = (double **)malloc(sizeof(double *) *subRows);

 /*for (int i = 0; i < subRows; i++) // for each process allocate rows of size N
  {
    u[i] = (double *)malloc(sizeof(double) * N);
    w[i] = (double *)malloc(sizeof(double) * N);
  }*/
 
  // double u[M][N];
  // double w[M][N];
  double wtime;

  // TODO: Make sure that only one MPI process prints the output (DONE)
  if (world_rank == 0)
  {
    std::cout << "\n";
    std::cout << "HEATED_PLATE\n";
    std::cout << "  C++ version\n";
    std::cout
        << "  A program to solve for the steady state temperature distribution\n";
    std::cout << "  over a rectangular plate.\n";
    std::cout << "\n";
    std::cout << "  Spatial grid of " << M << " by " << N << " points.\n";
    std::cout << "  The iteration will be repeated until the change is <= " << epsilon
         << "\n";
  
  // TODO: Get number of MPI processes (DONE)
  long num_procs = world_size; // num_procs = number of processes
  std::cout << "  Number of processes available = " << num_procs << "\n";
  }
  std::string refFile;

  if (argc >= 2)
  {
    refFile = argv[1];
    std::cout << "  Reading reference results from " << refFile << std::endl;
  }
  else
  {
    std::cout << "  No reference file provided - skipping verification." << std::endl;
  }

  //
  //  Set the boundary values, which don't change.
  // TODO: Be aware that each process hold a different part of the matrix and may only (ALMOST DONE, but problem with invalid numers of process)
  // initialize his own part
  mean = 0.0;

  // All process has the same subset number of rows to process

  if (world_rank == 0) // first process
  {
    //std::cout << "  Settings boundaries(master)" << std::endl;
    for (i = 0; i < N; i++) // first row has only zeros
    {
      w[0][i] = 0.0;
    }
    for (j = 1; j < subRows; j++) // for other rows of first process
    {
      w[j][0] = 100.0;     // left value is 100
      w[j][N - 1] = 100.0; // right value is 100
    }

  }

  if (world_rank < (world_size-1) && world_rank!=0) // all processes except the first and last one
  {
    //std::cout << "  Settings boundaries(middle worker)" << std::endl;
    for (j = 0; j < subRows; j++) // for j-th row, subRows = number of rows per process
    {
      w[j][0] = 100.0;     // left value is 100
      w[j][N - 1] = 100.0; // right value is 100
    }
    //std::cout << "  Set boundaries(middle worker)" << std::endl;
  }

  if (world_rank == world_size-1) // last process gives last row only 100.0
  {
    //std::cout << "  Settings boundaries(last process)" << std::endl;
    for (j = 0; j < subRows-1; j++) // for j-th row
    {
      w[j][0] = 100.0;     // left value is 100
      w[j][N - 1] = 100.0; // right value is 100
    }
    for (i = 0; i < N; i++) // last row has only 100
    {
      w[subRows-1][i] = 100.0;
    }
    //std::cout << "  Set boundaries(last process)" << std::endl;
  }

  //
  //  Average the boundary values, to come up with a reasonable
  //  initial value for the interior.
  //  TODO: Keep in mind that you need to average over all processes (DONE)

  for (i = 1; i < subRows - 1; i++)
  {
    mean = mean + w[i][0] + w[i][N - 1]; // sum of each row, except the first and the last one
  }
  for (j = 0; j < N; j++)
  {
    mean = mean + w[subRows - 1][j] + w[0][j]; // sum of first and last row
  }

  // Master process should receive all processes sums to calculate mean and send the same updated mean to all other processes
  double sum = 0.0;
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&mean, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (world_rank == 0)
  {   
    mean = sum / (double)(2 * M + 2 * N - 4);  
                  // mean value Temperature of matrix
    //MPI_Bcast(&mean, 1, MPI_DOUBLE, world_rank, MPI_COMM_WORLD); // Send to all other processes (world_rank)
  }
  MPI_Bcast(&mean, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(world_rank == 0)
  {
    std::cout << "\n";
    std::cout << "  MEAN = " << mean << "\n";
  }
  if(world_rank == 0){
    for(i=1;i<subRows;i++){
      for(j=1;j<N-1;j++){
        w[i][j] = mean;
      }
    }
  }
  if(world_rank < world_size-1 && world_rank != 0){
    for(i=0;i<subRows;i++){
      for(j=1;j<N-1;j++){
        w[i][j] = mean;
      }
    }
  } 
  if(world_rank == world_size-1){
    for(i=0;i<subRows-1;i++){
      for(j=1;j<N-1;j++){
        w[i][j] = mean;
      }
    }
  }
   //MPI_Barrier(MPI_COMM_WORLD);

  //
  //  Initialize the interior solution to the mean value.
  //  All matrix elements gets same temperature (mean)

  //
  //  Iterate until the  new solution W differs from the old solution U
  //  by no more than EPSILON.
  //
  iterations = 0;
  iterations_print = 1;
  if(world_rank==0){
    std::cout << "\n";
    std::cout << " Iteration  Change\n";
    std::cout << "\n";
    std::cout << "------------------" << std::endl;
  }


  Timer timer;
  double next_row[N];
  double prev_row[N];
  diff = epsilon;



  // TODO: This statement should yield the same result on all processes, so that all
  // processes stop at the same iteration
  while (epsilon <= diff)
  {
    //
    //  Save the old solution in U, when deviation exceeds epsilon
    //
    for (i = 0; i < subRows; i++)
    {
      for (j = 0; j < N; j++)
      {
        u[i][j] = w[i][j];
      }
    }

    //std::cout << "updated matrix" << std::endl;
    //
    //  Determine the new estimate of the solution at the interior points.
    //  The new solution W is the average of north, south, east and west
    //  neighbors.
    //  TODO: Here you may need parts of the matrix that are part of other processes
    if (world_rank == 0)    //Sends their last row to next process, receives last row from next process
    {
      MPI_Isend(&u[subRows-1][0], N, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD,&request);

      MPI_Irecv(&next_row[0], N, MPI_DOUBLE, world_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD,&request);


      for (i = 1; i < subRows -1 ; i++)
      {

        for (j = 1; j < N - 1; j++) //starting at first column, stops before last column
        {
          // w = (north + south + west + east)/4

          w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
        }
      }

      for (j = 1; j < N - 1; j++)
      {
        // w = (north + south + west + east)/4
        w[subRows-1][j] = (u[subRows - 2][j] + next_row[j] + u[subRows-1][j - 1] + u[subRows-1][j + 1]) / 4.0;
      }

    }

    if (world_rank != 0 && world_rank != world_size-1) //all other process except for first and last sends their first and last rows to the previous and next process respectively
    {                                                //receives first and last row from these processes as well

      MPI_Isend(&u[0][0],N,MPI_DOUBLE, world_rank-1,0,MPI_COMM_WORLD,&request);
      MPI_Isend(&u[subRows-1][0],N,MPI_DOUBLE, world_rank+1,0,MPI_COMM_WORLD,&request);

      MPI_Irecv(&prev_row[0],N,MPI_DOUBLE,world_rank-1,MPI_ANY_TAG,MPI_COMM_WORLD, &request);
      MPI_Irecv(&next_row[0],N,MPI_DOUBLE,world_rank+1,MPI_ANY_TAG,MPI_COMM_WORLD, &request);

      for (i = 1; i < subRows - 1; i++)
      {
        for (j = 1; j < N - 1; j++)
        {
          // w = (north + south + west + east)/4
          w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
        }
      }

      for (i = 0; i < N ; i++)
      {    
        // w = (north + south + west + east)/4
        w[0][i] = (prev_row[i] + u[1][i] + u[0][i - 1] + u[0][i + 1]) / 4.0;
        
      }

      for (i = 0; i < N ; i++)
      {    
        // w = (north + south + west + east)/4
        w[subRows-1][i] = (u[subRows-2][i] + next_row[i] + u[subRows-1][i - 1] + u[subRows-1][i + 1]) / 4.0;
        
      }

    }

    if (world_rank == world_size-1)   //last process sends first row to previous process and receives last row from previous process
    {

      MPI_Isend(&u[0][0], N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD,&request);
      MPI_Irecv(&prev_row[0], N, MPI_DOUBLE, world_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
      
      for (i = 0; i < N ; i++)     //first row
      {    
        // w = (north + south + west + east)/4
        w[0][i] = (prev_row[i] + u[1][i] + u[0][i - 1] + u[0][i + 1]) / 4.0;
        
      }

      for (i = 1; i < subRows - 1; i++) //all other rows except last
      {
        for (j = 1; j < N - 1; j++)
        {
          // w = (north + south + west + east)/4
          w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
        }
      }

    }
    diff = 0.0;
    partialDiff = 0.0;
    // TODO: Be aware that each process only computes its local diff here. You may
    // need to combine the results from all processes
    for (i = 1; i < subRows - 1; i++)
    {
      for (j = 1; j < N - 1; j++)
      {
        if (partialDiff < fabs(w[i][j] - u[i][j])) // if diff < 0, then update its current value
        {
          partialDiff = fabs(w[i][j] - u[i][j]);
        }
      }
    }
    MPI_Reduce(&partialDiff,&diff,1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&diff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(world_rank==0){
      iterations++;
      if (iterations == iterations_print)
      {
        std::cout << "  " << setw(8) << iterations << "  " << diff << "\n";
        iterations_print = 2 * iterations_print;
      }
    }

  }

  // TODO: Insert a Barrier before time measurement (DONE)
  
  MPI_Barrier(MPI_COMM_WORLD);
  wtime = timer.getSeconds();

  // TODO: Make sure that only one MPI process prints the output (DONE)
  if (world_rank == 0)
  {
    std::cout << "\n";
    std::cout << "  " << setw(8) << iterations << "  " << diff << "\n";
    std::cout << "\n";
    std::cout << "  Error tolerance achieved.\n";
    std::cout << "  Wallclock time = " << setprecision(3) << wtime << "s\n";
  
  //
  //  Terminate.
  //
  
    std::cout << "\n";
    std::cout << "HEATED_PLATE:\n";
    std::cout << "  Normal end of execution.\n";
  }
    if (!refFile.empty())
    {
      double ref[M][N];

        if (!readFromFile(refFile, &ref[0][0], M * N))
        {
          std::cout << "  Verification failed - could not load reference file." << std::endl;
          return 0;
        }
        if(world_rank==0){


          if (verify(&w[0][0], &ref[0][0], subRows*N))
          {
            std::cout << "  Result is valid." << world_rank <<std::endl;
          }
          else
          {
            std::cout << "  Verification failed!" << world_rank <<std::endl;
          }
        }
        else{

          if (verify(&w[0][0], &ref[(world_rank*subRows) + modulo][0],subRows*N))
          {
            std::cout << "  Result is valid." << world_rank <<std::endl;
          }
          else
          {
            std::cout << "  Verification failed!" << world_rank <<std::endl;
          }
        }

    
  }

  
  
  
  return 0;
#undef M
#undef N
MPI_Finalize();
}
