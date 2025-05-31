#include "utils.h"

void block_dist_vec(double* source_vec, double*& local_vec, MPI_Comm comm, int n){
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // row distribute
    MPI_Comm col_subcom;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims, &col_subcom);
    //at first col
    if(coords[1]==0){
        // calc count and displacement
        std::vector<int> send_counts(dims[0]);
        std::vector<int> send_displs(dims[0]);
        for(int i = 0; i<dims[0]; i++){
            send_counts[i] = find_block_num(n, dims[0], i);
        }
        send_displs[0] =0;
        for(int i = 1; i<dims[0]; i++){
            send_displs[i] = send_displs[i-1] + send_counts[i-1];
        }
        local_vec = new double[send_counts[coords[0]]];
        
        //get root of the col comm
        int rrank;
        int rcoords[1] = {0};
        MPI_Cart_rank(col_subcom, rcoords, &rrank);
        MPI_Scatterv(source_vec, &send_counts[0], &send_displs[0], MPI_DOUBLE, local_vec, send_counts[coords[0]], MPI_DOUBLE, rrank, col_subcom);
    }
    MPI_Comm_free(&col_subcom);
    return;
}

void block_dist_mat(double* source_mat, double*& local_mat, MPI_Comm comm, int n){
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    int local_mat_rows = find_block_num(n, dims[0], coords[0]);
    int local_mat_cols = find_block_num(n, dims[1], coords[1]);
    local_mat = new double[local_mat_rows*local_mat_cols];
    
    // row distribute
    MPI_Comm col_subcom;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims, &col_subcom);
    double* row_mat = NULL; // matrix block after row distribute
    // at first col
    if(coords[1]==0){
        // calc count and displacement
        std::vector<int> send_counts(dims[0]);
        std::vector<int> send_displs(dims[0]);
        for(int i = 0; i<dims[0]; i++){
            // note n*(n/p) data per row block
            send_counts[i] = find_block_num(n, dims[0], i) * n;
        }
        send_displs[0] =0;
        for(int i = 1; i<dims[0]; i++){
            send_displs[i] = send_displs[i-1] + send_counts[i-1];
        }
        row_mat = new double[send_counts[coords[0]]];
        //get root of the col comm
        int rrank;
        int rcoords[1] = {0};
        MPI_Cart_rank(col_subcom, rcoords, &rrank);
        MPI_Scatterv(source_mat, &send_counts[0], &send_displs[0], MPI_DOUBLE, row_mat, send_counts[coords[0]], MPI_DOUBLE, rrank, col_subcom);
    }

    // col distribute
    MPI_Comm row_subcom;
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm, remain_dims, &row_subcom);
    
    // calc count and displacement
    std::vector<int> send_counts2(dims[1]);
    std::vector<int> send_displs2(dims[1]);
    for(int i = 0; i<dims[1]; i++){
            send_counts2[i] = find_block_num(n, dims[1], i);
    }
    send_displs2[0] =0;
    for(int i = 1; i<dims[1]; i++){
        send_displs2[i] = send_displs2[i-1] + send_counts2[i-1];
    }
    
    //get root of the row comm
    int rrank;
    int rcoords[1] = {0};
    MPI_Cart_rank(row_subcom, rcoords, &rrank);
    for(int i=0; i<local_mat_rows; i++){
        MPI_Scatterv(row_mat+i*n, &send_counts2[0], &send_displs2[0], MPI_DOUBLE, local_mat+i*local_mat_cols, local_mat_cols, MPI_DOUBLE, rrank, row_subcom);
    }

    delete [] row_mat;
    MPI_Comm_free(&col_subcom);
    MPI_Comm_free(&row_subcom);
    return;
}

void diagnal_bcast_vec(double* source_vec, double* dest_vec, MPI_Comm comm, int n){

    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    int row_block_num = find_block_num(n, dims[0], coords[0]);
    int col_block_num = find_block_num(n, dims[1], coords[1]);

    // col and row subcommunicators
    MPI_Comm col_comm, row_comm;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims, &col_comm);
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    // sending and receiving
    // first row sending to diagonal
    if (coords[1] == 0 && coords[0] != 0){
        int dest_coords[2] = {coords[0], coords[0]};
        // find rank of this dest_coords
        int dest_rank;
        MPI_Cart_rank(comm, dest_coords, &dest_rank);
        MPI_Send(source_vec, row_block_num, MPI_DOUBLE, dest_rank, 0, comm);
    }
    // diagonal receiving 
    if (coords[1] == coords[0] && coords[0] != 0){
        int sour_coords[2] = {coords[0], 0};
        // find rank of this sour_coords
        int sour_rank;
        MPI_Cart_rank(comm, sour_coords, &sour_rank);
        MPI_Recv(dest_vec, row_block_num, MPI_DOUBLE, sour_rank, 0, comm, MPI_STATUS_IGNORE);
    }
    // for the left top block
    if (coords[1] == 0 && coords[0] == 0){
        for (int k = 0; k < col_block_num; k++){
            dest_vec[k] = source_vec[k];
        }
    }
    
    // broadcast values on diagonal to the whole grid
    // done by broadcast value on (j, j) to col j
    MPI_Bcast(dest_vec, col_block_num, MPI_DOUBLE, coords[1], col_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);

}

void p_matrix_vector_multi(double* local_A, double* local_x, double* local_y, MPI_Comm comm, int n){
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    //get number of elements in row and col
    int local_mat_rows = find_block_num(n, dims[0], coords[0]);
    int local_mat_cols = find_block_num(n, dims[1], coords[1]);

    //bcast local x to all process
    double* local_bcasted_x = new double[local_mat_cols];
    diagnal_bcast_vec(local_x, local_bcasted_x, comm, n);

    //calc local temp y value before reduce
    double* local_temp_y = new double[local_mat_rows];
    matrix_vector_multi(local_mat_rows, local_mat_cols, local_A, local_bcasted_x, local_temp_y);

    //find the row comm and reduce to rank 0
    MPI_Comm row_subcom;
    int remain_dims[2] = {0, 1};
    MPI_Cart_sub(comm,remain_dims, &row_subcom);
    int rrank;
    int rcoords[1] = {0};
    MPI_Cart_rank(row_subcom, rcoords, &rrank);
    MPI_Reduce(local_temp_y, local_y, local_mat_rows, MPI_DOUBLE, MPI_SUM, rrank, row_subcom);

    delete [] local_bcasted_x;
    delete [] local_temp_y;
    MPI_Comm_free(&row_subcom);
    return;
}

void p_jacobi(double* local_A, double* local_b, double* local_x, MPI_Comm comm, int n, int max_iteration, double epsilon){
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    int local_mat_rows = find_block_num(n, dims[0], coords[0]);
    int local_mat_cols = find_block_num(n, dims[1], coords[1]);

    // col and row subcommunicators
    MPI_Comm col_subcom, row_subcom;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims, &col_subcom);
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm, remain_dims, &row_subcom);

    // Initialize x
    if (coords[1] == 0){
        for (int i = 0; i < local_mat_rows; i++) {
            local_x[i] = 0;
        }
    }


    // Initialize R
    double* local_R = new double[local_mat_rows*local_mat_cols];
    std::copy(local_A, local_A+local_mat_rows*local_mat_cols, local_R);
    if (coords[0] == coords[1]){
        for (int i = 0; i < local_mat_rows; i++){
            local_R[i+i*local_mat_cols] = 0;
        }
    }

    // Initialize D
    std::vector<double> local_D(local_mat_cols);
    // diagnal processor calc local_D value
    if (coords[0] == coords[1]){
        for (int i = 0; i < local_mat_rows; i++){
            local_D[i] = local_A[i+i*local_mat_cols];
        }
        // send to first col 
        // exclude (0,0), no need to send
        if(coords[0]!=0){
            int rrank;
            int rcoords[1] = {0};
            MPI_Cart_rank(row_subcom, rcoords, &rrank);
            MPI_Send(&local_D[0],local_mat_rows,MPI_DOUBLE, rrank, 1, row_subcom);
        } 
    }
    
    //receive at first col
    // but (0,0) does not need to receive anything, its local_D has value
    if(coords[0]!=0 && coords[1]==0){
        int sender_rank;
        int sender_coords[1] = {coords[0]};// in row_subcom, diagnal of row i is at ith index
        MPI_Status stat;
        MPI_Cart_rank(row_subcom, sender_coords, &sender_rank);
        MPI_Recv(&local_D[0], local_mat_rows, MPI_DOUBLE, sender_rank, 1, row_subcom, &stat);
    }

    for(int i =0; i<max_iteration; i++){
        // Compute Rx 
        std::vector<double> local_temp(local_mat_rows);
        p_matrix_vector_multi(local_R,local_x,&local_temp[0],comm,n);

        // Compute b - Rx
        // Compute x =  D^{-1}(b-Rx)
        if (coords[1] == 0){
            for (int j = 0; j < local_mat_rows; j++){
                local_temp[j] = local_b[j] - local_temp[j];
                local_x[j] = local_temp[j]/local_D[j];
            }
        }

        // Compute Ax
        std::vector<double> local_y(local_mat_rows);
        p_matrix_vector_multi(local_A,local_x,&local_y[0],comm,n);

        // Compute Ax-b
        // Compute ||Ax-b||
        int terminate = 0;
        double norm_sum = 0.0;
        double norm = 0.0;
        if (coords[1] == 0){
            for (int j = 0; j < local_mat_rows; j++){
                norm+=std::pow((local_y[j] - local_b[j]), 2);
            }
            MPI_Allreduce(&norm, &norm_sum, 1, MPI_DOUBLE, MPI_SUM, col_subcom);
            if (sqrt(norm_sum) <= epsilon){
                terminate = 1;
            }
        }
        int rrank;
        int rcoords[1] = {0};
        MPI_Cart_rank(row_subcom, rcoords, &rrank);
        MPI_Bcast(&terminate, 1, MPI_INT, rrank, row_subcom);
        if (terminate == 1){
            delete[] local_R;
            return;
        }
    }
    delete[] local_R;
    MPI_Comm_free(&col_subcom);
    MPI_Comm_free(&row_subcom);
    return;
}

void gather_vec(double* source_vec, double* dest_vec, MPI_Comm comm, int n){
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    //col communicator
    MPI_Comm col_subcom;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims, &col_subcom);
    
    //in first col
    if(coords[1]==0){
        std::vector<int> rec_counts(dims[0]);
        std::vector<int> rec_displs(dims[0]);
        for(int i=0; i<dims[0];i++){
            rec_counts[i] = find_block_num(n, dims[0], i);
        }
        rec_displs[0] = 0;
        for(int i=1; i<dims[0]; i++){
            rec_displs[i] = rec_displs[i-1]+rec_counts[i-1];
        }

        //get root rank
        int rrank;
        int rcoords[1] = {0};
        MPI_Cart_rank(col_subcom, rcoords, &rrank);
        MPI_Gatherv(source_vec, rec_counts[coords[0]], MPI_DOUBLE, dest_vec, &rec_counts[0], &rec_displs[0], MPI_DOUBLE, rrank, col_subcom);
    }
    MPI_Comm_free(&col_subcom);
    return;
}

int main(int argc, char* argv[]){
    // Intialize and handle errors
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    int q = (int)sqrt(world_size);
    
    MPI_Comm grid;
    int dims[2] = {q,q};
    int periods[2] = {0,0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims , periods, 1, &grid);

    // Read in n, A, b.
    int n;
    double* A = NULL;
    std::vector<double> b;

    int crank0;
    int rcoords[2] = {0,0};
    MPI_Cart_rank(grid, rcoords, &crank0);
    if (world_rank == crank0){
        read_matrix(argv[1],n, A);
        b = read_vector(argv[2],n);
    }

    
    // Let all processors know n.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Block distribute A and b
    double* local_A= NULL;
    double* local_b =NULL;
    block_dist_mat(A, local_A, grid, n);
    block_dist_vec(&b[0], local_b, grid, n);

    // Parallel Jacobi, runtime
    double time_init,time_end;
    time_init = MPI_Wtime();
    int local_x_size = find_block_num(n, q, 0);
    std::vector<double> local_x(local_x_size);
    int max_iter = 1000000;
    double epsilon = 1e-9;
    p_jacobi(local_A, local_b, &local_x[0], grid, n, max_iter, epsilon);
    time_end = (MPI_Wtime() - time_init)*1000;
    
    // Gather x vector
    std::vector<double> x(n);
    gather_vec(&local_x[0],&x[0],grid,n);

    // Write x vector and run time
    if (world_rank == crank0){
        write_vector(argv[3],x,false,time_end);
    }
    MPI_Comm_free(&grid);
    MPI_Finalize();
    return 0;
}
