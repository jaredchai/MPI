#include <mpi.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

int find_block_num(int n, int p, int rank){
    return n/p +((rank < n%p)? 1:0);
}

void read_matrix(const std::string filename, int& n, double*& matrix){
    std::ifstream input_file(filename);

    if (!input_file) {
        std::cerr << "Error: Unable to open the file" << std::endl;
        return;
    }

    input_file >> n; // Read the value of n
    matrix = new double[n*n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            input_file >> matrix[i*n+j]; // Read matrix elements
        }
    }

    input_file.close();

    return;
}

std::vector<double> read_vector(const std::string filename, const int n){
    std::ifstream input_file(filename);
    std::vector<double> vector(n);

    if (!input_file) {
        std::cerr << "Error: Unable to open the file" << std::endl;
        return vector;
    }


    for (int i = 0; i < n; i++) {
        input_file >> vector[i];
    }

    input_file.close();

    return vector;
}

void write_vector(const std::string filename, std::vector<double> x,
    bool add_time = false, double timing = 0.0){
    std::ofstream output_file(filename);

    if (!output_file) {
        std::cerr << "Error: Unable to open the file" << std::endl;
        return;
    }

    if(add_time){
        output_file << std::fixed << std::setprecision(6) << timing << std::endl;
    }


    for (int i=0; i<x.size(); i++) {
        output_file<< std::fixed << std::setprecision(16) << x[i] << " "; // Write vector elements
    }

    output_file.close();
}

// mxn A
// nx1 x
// mx1 y
void matrix_vector_multi(const int m, const int n, const double* A, const double* x, double* y){
    for (int i = 0; i<m; i++){
        y[i] = 0.0;
        for (int j = 0; j<n; j++){
            y[i] += A[i*n + j]*x[j];
        }
    }
}