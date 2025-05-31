Run the program with `qsub int_calc.pbs`. Or run `mpirun -np [process number] ./int_calc  [argument n]`

The cpp code is in `int_calc.cpp`. Compile the exe `int_calc` by `make`.

The code has two helper function `f(x)` that return the function value and `int_app(lower, upper, h)` that return the approximation of the integration from `lower` to `upper`. The `main` functio will first calculate its local integration the use `MPI_Reduce` to sum all local value to processer 0.

`int_calc.out` will contain the output of the code.
