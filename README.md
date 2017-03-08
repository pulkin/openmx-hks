# openmx-hks

Parse the tight-binding Hamiltonian and other data from [OpenMX](http://www.openmx-square.org/)
and transform/output it into various formats.

## Features

1. Portable: only `*.hks` file is required
2. Easy: the tool is not compiled across huge OpenMX libraries
3. Functional:
   - all data is parsed including the Hamiltonian, the overlap matrix,
     the Hartree potential, the Fermi level, the crystal structure and more;
   - capable of changing and copying the Fermi level;
   - shift the Hamiltonian by a constant;
   - save the Hamiltonian and the overlap matrices into a mat or json file;
  
## Downloading and compiling

No special libraries needed, following should work:
```
mkdir openmx-hks
git clone https://github.com/pulkin/openmx-hks.git
cd openmx-hks/src
make
```

## Running

1. Run your OpenMX calculation with options `NEGF.output_hks on` and
   `NEGF.filename.hks  my_hks_file.hks`
2. Run `openmx-hks` to display information `./openmx-hks display path_to_my_hks_file.hks`
   or create a Matlab-formatted file `./openmx-hks extract-hamiltonian path_to_my_hks_file.hks matlab_file.mat`
3. Use this data in Matlab (or everywhere else) to plot band structures,
   create effective Hamiltonians, study other properties not accessable
   directly from OpenMX.

## Bugs

Report [here](https://github.com/pulkin/openmx-hks/issues).
