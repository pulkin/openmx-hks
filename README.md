[![Build Status](https://dev.azure.com/gpulkin/openmx-hks/_apis/build/status/pulkin.openmx-hks?branchName=master)](https://dev.azure.com/gpulkin/openmx-hks/_build/latest?definitionId=1&branchName=master)

# openmx-hks

Parse the tight-binding Hamiltonian and other data from [OpenMX](http://www.openmx-square.org/)
and transform/output it into various formats.
Examples also include plotting
electronic band structure and calculating ballistic transport properties.

## Features

1. Portable: only `*.hks` file is required
2. Easy: minimal dependencies, straightforward options
3. Functional:
   - all data is parsed including the Hamiltonian, the overlap matrix,
     the Hartree potential, the Fermi level, the atomic structure and more;
   - adjusts energies;
   - exports Hamiltonian and overlap matrices into various formats
     (MATLAB, json, hdf5);
   - exports structure into an XSF file;
  
## Download

From the [releases page](https://github.com/pulkin/openmx-hks/releases/tag/latest-build)

## Compiling

Ubuntu example

Install dependencies
```
sudo apt-get install build-essential gcc-multilib libhdf5-dev
sudo ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial.so /usr/lib/x86_64-linux-gnu/libhdf5.so
```
Clone and make
```
mkdir openmx-hks
git clone https://github.com/pulkin/openmx-hks.git
cd openmx-hks
make -C src
```

## Running examples

```
cd openmx-hks/examples/01_hks2json
./run
```

## Bugs

Report [here](https://github.com/pulkin/openmx-hks/issues)
