# Metropolis algorithm for the 3D-ising ising model

## Code Deployment

### Compilation

To compile the c++ code you can use the standard cmake-make build. To build with code optimization in this folder the commands are

```
cmake -DCMAKE_BUILD_TYPE=Release .
make -j
```

### Execution

The compilation will produce a first executable called `ising`. You can run it by passing an input file as an argument. The default is `input.json`.
The executable will produce in the folder `outputs` a file with the magnetization and energies sample for each value of beta and L requested in the input file. The Wolff algoeithm is used.

The second executable is called `autocorrelation` and will compute the execution time and autocorrelation for hardcoded temperatures and cluster sizes, for both the Metropolis and the Wolff algorithm.

### Analysis

The provided python scripts plots the Binder cumulants. No arguments is required and they assume the data is the `output` folder and
`input.json` contains the desired cluster sizes.
Feel free to modify them if different quantities are to be measured or you want to use data in a different folder.

## Discussion
From an analysis of the binder cumulant relative to the magnetization, we can validate the correctness of the Wolff algorithm and obtain. in just a minute of runtime on an Intel i7-8700K, the data
relative to the included figure `binder_cumulant.pdf`. We can observe that the estimate of the critical temperature is more accurate than what we obtained with the Metropolis code.

We can compare the two codes also by looking at the run-time and autocorrelation time, reported in the tables `timing_<algorithm>.txt`. We can see that there is not a large difference for high temperature, where the clusters are small and Wolff resembles a single spin update, albeit always accepted. We can observe the advantage of the cluster algorithm more clearly near criticality or in the ferromagnetic phase. The autocorrelation time of the Wolff algorithm at critical temperature still scales as a powe law with the lattice size, but the critical exponent compared to the Metropolis algorithm is smaller. (Quantitative analysis omitted here).

### Implementation details.
The Wolff and Metropolis algorithm use the C++ inherithance feature to share the common interface among them. Storing the each lattice site in a single byte integer proved to be beneficial compared to
a bitmap (simply reling on `std::vector<bool>) due to the associated bit manipulation cost. Mind that more refined techniques where an ensamble of 8 lattices is represented in a single byte and updated in parallel might be faster, but more cumberome to implement.

The Wolff algorithm performs a breadth first search of the cluster and flips it on the fly. Mind that it is not necessary to keep a list of visited sites, as fliped spins can not reenter in the cluster search.

