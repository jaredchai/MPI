#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<mpi.h>

// Fucntion f(x) to be integrated
double f(double x){
    return 4/(1+x*x);
}

// Integral approximation from lower to upper with h per step
double int_app(double lower, double upper, double h){
    double int_sum;
    double x;
    for (x = lower; x< upper; x+=h){
        int_sum += f(x+0.5*h)*h;
        // printf("integrate sum : %f, input: %f\n", int_sum, x-0.5*h);
    }
    return int_sum;
}

int main(int argc, char *argv[]){
    // set up argument
    int n = atoi(argv[1]);
    double h = 1.0/n;
    double start, end;
    
    //initailize MPI
    MPI_Init(&argc,&argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    

    //p not divide by n
    int remainder = n%world_size;
    if(remainder != 0){
        n -= remainder;
        h = 1.0/n;
    }

    //start timer
    if(world_rank == 0){
        start = MPI_Wtime();
    }

    //integrate locally
    double local_int_sum = 0;
    double local_interval = n/world_size*h;
    double local_lower = world_rank*local_interval;
    double local_upper = local_lower+local_interval;
    //printf("rank %d: integrate from %f to %f with delta %f\n", world_rank, local_lower, local_upper, h);
    local_int_sum = int_app(local_lower, local_upper, h);

    //reduce the global intergral sum to root
    double global_int_sum;
    MPI_Reduce(&local_int_sum, &global_int_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(world_rank == 0){
        end = MPI_Wtime();
        printf("%.11f, %f\n", global_int_sum, end-start);
    }

    MPI_Finalize();
    return 0;

}