# Installation

- The NVIDIA CUDA compiler needs to be installed: check if it is so by running `which nvcc` in your shell. If nothing comes up, follow the instructions on the NVIDIA developer page at
https://developer.nvidia.com/cuda-downloads. Have fun with that one!!

- This code has an extra dependency for reading the input file, `PTL`. The repository is located at https://github.com/wvannoordt/PTL.git.
To build this, simple clone the repository somewhere and type `make`. Afterwards, the environment variable `$PTL` needs to be defined and
point to the main directory in the repo, e.g. put `export PTL='/path/to/PTL'` in your bashrc. type `echo $PTL` to make sure you have it defined,
you may need to run `source ~/.bashrc`.

- Once `PTL` has been installed and the CUDA compilers are installed, simple type `make` in the main directory. This will create an executable in the `exampleRun` directory called
`gpu-tgv`. This can now be run without any arguments.

- Visualization of the data can be done using `VisIt`, which can be downloaded directly from their website at https://wci.llnl.gov/simulation/computer-codes/visit/executables