import numpy as np
import sys
import os

def generate_A(n):
    a = np.random.randn(n, n)
    a = a + np.linalg.norm(a, 1)*np.identity(n)
    return a


def generate_x(n):
    x = np.random.randn(n)
    return x


def write_to_file(n, name):
    A = generate_A(n)

    # for band matrix
    # for i in range(n):
    #     for j in range(n):
    #         if(i-j>1 or i-j < -1):
    #             A[i][j] = 0

    x = generate_x(n)
    b = A @ x
    A = A.flatten()
    if not os.path.exists("./input"):
        os.makedirs("./input")
    if not os.path.exists("./output"):
        os.makedirs("./output")
    with open("./input/input_mat_" + str(name) + ".txt", 'w') as f:
        f.write(str(n) + "\n")
        for j in range(n**2):
            f.write(str(A[j]) + " ")
            if j % n == n-1:
                f.write("\n")
    f.close()
    with open("./input/input_vec_" + str(name) + ".txt", 'w') as f:
        for j in range(n):
            f.write(str(b[j]) + " ")
    f.close()
    with open("./output/out_vec_" + str(name) + ".txt", 'w') as f:
        for j in range(n):
            f.write(str(x[j]) + " ")
    f.close()


def main():
    # if len(sys.argv) < 2:
    #     print("Usage: python generate_input.py <n>")
    #     return
    n = int(sys.argv[1])
    name = sys.argv[2]
    write_to_file(n, name)


if __name__ == '__main__':
    main()