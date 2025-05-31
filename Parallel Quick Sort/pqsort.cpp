#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<algorithm>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<numeric>
#include<iomanip>
#include<random>

void qsort_local(std::vector<int> &local_v);
int find_pivot(int seed, std::vector<int> &local_v, int q, int m, MPI_Comm comm);
int partition(std::vector<int> &local_v, int pivot); 
int redistribute_data(std::vector<int> &local_v, int q, int m, int left_size, int right_size, MPI_Comm comm, std::vector<std::vector<int> > &all_to_all_send_counts, int &color);
void qsort_para(std::vector<int> &local_v, int size, MPI_Comm comm);

// serial sort the local vector when q = 1
void qsort_local(std::vector<int> &local_v){
    std::sort(local_v.begin(), local_v.end());
}

// use the seed to generate a random number from 0 to m-1
// return the pivot value at the kth number
int find_pivot(int seed, std::vector<int> &local_v, int q, int m, MPI_Comm comm){
    // seed the generator
    std::mt19937 gen(seed);
    // generate k between 0 and m-1;
    std::uniform_int_distribution<> dis(0, m - 1);
    int k = dis(gen);
    // compute the minimum number of elements in each process and the number of processes with the maximum number of elements
    int minSize =  m / q;
    int numMore = m % q;
    int maxSize;
    if (numMore == 0){
        numMore = q;
        maxSize = minSize;
    }
    else{
        maxSize = minSize+1;
    }
    // compute the community rank of the processor that holds the k-th element, as well as the index of this k-th element within this procesor.
    int targetRank, index;
    if ((numMore * maxSize) >= k){
        targetRank = k / maxSize;
        index = k - targetRank * maxSize;
    }
    else {
        targetRank = numMore + (k - (numMore * maxSize)) / minSize;
        index = k - (numMore * maxSize) - (targetRank - numMore) * minSize;
    }
    // broacast the value of the k-th element to all processors from the processor that holds it.
    int pivotValue, commRank;
    MPI_Comm_rank(comm, &commRank);
    if (commRank == targetRank){
        pivotValue = local_v[index];
    }
    MPI_Bcast(&pivotValue, 1, MPI_INT, targetRank, comm);
    return pivotValue;
}
// use the pivot value to partition the local vector for each processor
// return the index of the partition line <=> number of values less or equal pivot
int partition(std::vector<int> &local_v, int pivot){
    int i = 0, j = local_v.size() - 1;
    while (i <= j) {
        while (i <= j && local_v[i] <= pivot) i++;
        while (i <= j && local_v[j] > pivot) j--;
        if (i <= j) {
            int temp = local_v[i];
            local_v[i] = local_v[j];
            local_v[j] = temp;
            i++;
            j--;
        }
    }
    return i;
}

// the most challenging part
// 1. All Gather partition numbers m' and m''
// 2. find the new processor ratios of lesser and greater group, round to int
// 3. within each group, cylinder distribute data for load balance
// 4. calculate send_counts, send_displs, receive_counts, receive_displs 
// 5. AlltoAllv
// 6. pass a group id for comm split
int redistribute_data(std::vector<int> &local_v, int q, int m, int left_size, int right_size, MPI_Comm comm, std::vector<std::vector<int> > &all_to_all_send_counts, int &color){
    // All gather into an array of size 2 * number of processors in community
    // The array has format of [left size in proc 0; right size in proc 0; left size in proc 1; right size in proc 1;...]
    int partition[2] = {left_size,right_size};
    int partitionResults[2*q];
    MPI_Allgather(&partition, 2, MPI_INT, &partitionResults, 2, MPI_INT, comm);
    // Compute the sum of elements in the left and right of the procs. The two ints sum to m;
    int left_total = 0,right_total = 0;
    for (int i = 0; i < q; i++){
        left_total+=partitionResults[0+i*2];
    }
    right_total = m - left_total;
    // Compute the ratio "left_proc : right proc" so that it is close to "left_total : right_total" and sum to q (number of procs in comm);
    double r = (double)left_total / ((double)left_total+(double)right_total);
    double R = (double)left_total / (double)right_total;
    int round_down = std::floor(r*q);
    int round_up = std::ceil(r*q);
    double R1 = double(round_down)/double(q-round_down);
    double R2 = double(round_up)/double(q-round_up);
    int left_procs;
    if (std::abs(R1-R) > std::abs(R2-R)){
        left_procs = round_up;
    }
    else{
        left_procs = round_down;
    }
    int right_procs = q - left_procs;
    // Deal with the situation with 0 proc on one side.
    if (left_procs == 0 && right_procs > 0){
        left_procs+=1;
        right_procs-=1;
    }
    else if (left_procs > 0 && right_procs == 0){
        left_procs-=1;
        right_procs+=1;
    }
    else if (left_procs == 0 && right_procs == 0){
        printf("error. shouldn't happen.");
        return -1;
    }
    // fill the all to all send count table with cylinder distribution method
    int leftproc_ptr = 0, rightproc_ptr = left_procs;
    for (int i = 0; i < q; i++){
        // number to send to left side
        if (partitionResults[0+i*2] != 0){
            int minSize =  partitionResults[0+i*2] / left_procs;
            int numMore = partitionResults[0+i*2] % left_procs;
            int maxSize, numLess = 0;
            if (numMore == 0){
                numMore = left_procs;
                maxSize = minSize;
            }
            else{
                maxSize = minSize+1;
            }
            if (minSize != maxSize && minSize != 0){
                numLess = int((partitionResults[0+i*2] - numMore * maxSize)/minSize);
            }
            for (int j = 0; j < numMore;j++){
                all_to_all_send_counts[i][leftproc_ptr] = maxSize;
                leftproc_ptr++;
                if (leftproc_ptr >= left_procs){
                    leftproc_ptr = 0;
                }
            }
            int leftproc_subptr = leftproc_ptr;
            for (int j = 0; j < numLess;j++){
                all_to_all_send_counts[i][leftproc_subptr] = minSize;
                leftproc_subptr++;
                if (leftproc_subptr >= left_procs){
                    leftproc_subptr = 0;
                }
            }
        }
        // number to send to right side
        if (partitionResults[1+i*2] != 0){
            int minSize =  partitionResults[1+i*2] / right_procs;
            int numMore = partitionResults[1+i*2] % right_procs;
            int maxSize, numLess = 0;
            if (numMore == 0){
                numMore = right_procs;
                maxSize = minSize;
            }
            else{
                maxSize = minSize+1;
            }
            if (minSize != maxSize && minSize != 0){
                numLess = int((partitionResults[1+i*2] - numMore * maxSize)/minSize);
            }
            for (int j = 0; j < numMore;j++){
                all_to_all_send_counts[i][rightproc_ptr] = maxSize;
                rightproc_ptr++;
                if (rightproc_ptr >= q){
                    rightproc_ptr = left_procs;
                }
            }
            int rightproc_subptr = rightproc_ptr;
            for (int j = 0; j < numLess;j++){
                all_to_all_send_counts[i][rightproc_subptr] = minSize;
                rightproc_subptr++;
                if (rightproc_subptr >= q){
                    rightproc_subptr = left_procs;
                }
            }
        }   
    }
    // Compute the four vectors in Alltoallv;
    int commRank;
    MPI_Comm_rank(comm, &commRank);
    std::vector<int> send_counts = all_to_all_send_counts[commRank],receive_counts(q),receive_displs(q),send_displs(q);
    receive_counts[0] = all_to_all_send_counts[0][commRank];
    for (int i = 1; i < q; i++){
        receive_counts[i] = all_to_all_send_counts[i][commRank];
        receive_displs[i] = receive_displs[i-1]+all_to_all_send_counts[i-1][commRank];
        send_displs[i] = send_displs[i-1] + all_to_all_send_counts[commRank][i-1];
    }
    // Alltoallv
    std::vector<int> new_v(std::accumulate(receive_counts.begin(), receive_counts.end(), 0));
    MPI_Alltoallv(&local_v[0],&send_counts[0],&send_displs[0],MPI_INT,&new_v[0],&receive_counts[0],&receive_displs[0],MPI_INT,comm);
    local_v.resize(new_v.size());
    local_v = new_v;
    // Compute the color based on split of procs.
    color = (commRank >= left_procs);
    // Return the number of elements in the new community.
    if (color){
        return right_total;
    }
    else{
        return left_total;
    }
}   

// recursive pqsort
// 1. find pivot
// 2. partition data
// 3. redistribute data to left/right group processors
// 4. split comm and call recursive pqsort again
// base case: q=1, qsort_local()
void qsort_para(std::vector<int> &local_v, int m, MPI_Comm comm){
    // get current community q
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1){
        qsort_local(local_v);
        return;
    }
    MPI_Comm_rank(comm, &comm_rank);
    int newm = -1;
    int seed = 1;
    int color;
    // if one side is empty (selected min or max), redo the partition with increased seed;
    std::vector<int> local_v_cpy = local_v;
    int left_size,pivot;
    while (newm == -1){
        // find value of the pivot
        local_v = local_v_cpy;
        pivot = find_pivot(seed, local_v, comm_size, m, comm);
        // partition, find final index of the pivot;
        left_size = partition(local_v, pivot);
        // redistribute data
        std::vector<std::vector<int> > all_to_all_send_counts(comm_size,std::vector<int>(comm_size));
        newm = redistribute_data(local_v, comm_size, m, left_size, (local_v.size() - left_size), comm, all_to_all_send_counts, color);
        seed++;
    }
    MPI_Comm newcomm;
    MPI_Comm_split(comm, color, comm_rank, &newcomm);
    // Compute number of processor in the new community
    qsort_para(local_v, newm, newcomm);
}


//1. read input data to rank 0
//2. Bcast to all other processors, timer start
//3. call qsort_para()
//4. Gatherv all sorte data, timer end
//5. write output file
int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    int m;
    std::vector<int> local_v;
    std::vector<int> all_data;
    // read-in input file on proc 0
    if (world_rank == 0){
        std::ifstream infile(argv[1]);
        std::string line;
        std::getline(infile, line);
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        m = std::stoi(line);
        all_data.reserve(m);
        std::getline(infile, line);
        infile.close();
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        std::istringstream iss(line);
        std::string(num);
        while (iss >> num){
            all_data.push_back(std::stoi(num));
        }
    }
    // broadcast the total number of data to all procs
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int color = 0;
    MPI_Comm new_world;
    if (m < world_size && world_rank >= m ){
        color = 1;
    }
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &new_world);
    MPI_Comm_rank(new_world,&world_rank);
    MPI_Comm_size(new_world,&world_size);
    if (color == 0){   
        // compute the distribution of data
        int minSize =  m / world_size;
        int maxSize;
        int numMore = m % world_size;
        if (numMore == 0){
            numMore = world_size;
            maxSize = minSize;
        }
        else{
            maxSize = minSize+1;
        }
        std::vector<int> send_counts(world_size),send_disps(world_size);
        send_disps[0] = 0;
        std::fill(send_counts.begin(), (send_counts.begin() + numMore),maxSize);
        std::fill((send_counts.begin() + numMore), send_counts.end(), minSize);
        // scatter the data to the procs according to the distribution
        for (int i = 1; i < world_size; i++){
            send_disps[i] = send_disps[i-1]+send_counts[i-1];
        }
        int receiveCounts = send_counts[world_rank];
        local_v.resize(receiveCounts);
        MPI_Scatterv(&all_data[0],&send_counts[0],&send_disps[0],MPI_INT,&local_v[0],receiveCounts,MPI_INT,0,new_world);
        // start timer
        double time_init,time_end;
        time_init = MPI_Wtime();
        // Parallel quicksort
        qsort_para(local_v, m, new_world);

        // end timer
        time_end = (MPI_Wtime() - time_init)*1000;
        // Gather number of data on each proc
        std::vector<int> receive_counts(world_size), receive_disps(world_size);
        int sendData = local_v.size();
        MPI_Gather(&sendData,1,MPI_INT,&receive_counts[0],1,MPI_INT,0,new_world);
        // Gather the data from each proc to 0
        if (world_rank == 0){
            receive_disps[0] = 0;
            for (int i = 1; i < world_size; i++){
                receive_disps[i] = receive_disps[i-1]+receive_counts[i-1];
            }
        }
        int receiveSize = local_v.size();
        MPI_Gatherv(&local_v[0],receiveSize,MPI_INT,&all_data[0],&receive_counts[0],&receive_disps[0],MPI_INT,0,new_world);
        // Write to output file on proc 0
        if (world_rank == 0){
            std::ofstream outfile(argv[2]);
            for (int i = 0; i < m; i++){
                outfile << all_data[i] << " ";
            }
            outfile << "\n";
            outfile << std::setprecision(6) << std::fixed << time_end <<"\n";
            outfile.close();
        }
    }
    MPI_Finalize();
    return 0;
}
