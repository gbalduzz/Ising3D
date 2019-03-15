# Metropolis algorithm for the 3D-ising ising model

## Compilation

To compile the c++ code you can use the standard cmake-make build. To build with code optimization in this folder the commands are

```
cmake -DCMAKE_BUILD_TYPE=Release .
make -j
```

## Execution

The compilation will produce an executable called `ising`. You can run it by passing an input file as an argument. The default us `input.json`.
The executable will produce in the folder `outputs` a file with the magnetization and energies sample for each value of beta and L requested in the input file.

## Analysis

The provided python scripts plots the Binder cumulants and the magnetic susceptibiity. No arguments is required and they assume the data is the `output` folder and
`input.json` contains the desired cluster sizes.
Feel free to modify them if different quantities are to be measured or you want to use data in a different folder.