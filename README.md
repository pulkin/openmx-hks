# openmx-hks

Parse the tight-binding Hamiltonian and other data from [OpenMX](http://www.openmx-square.org/)
and transform/output it into various formats. Examples also include plotting
electronic band structure and calculating ballistic transport properties.

## Features

1. Portable: only `*.hks` file is required
2. Easy: the tool is not compiled across huge OpenMX libraries
3. Functional:
   - all data is parsed including the Hamiltonian, the overlap matrix,
     the Hartree potential, the Fermi level, the atomic structure and more;
   - adjusts energies;
   - exports Hamiltonian and overlap matrices into various formats
     (MATLAB, json);
   - exports structure into an XSF file;
  
## Downloading and compiling

No special libraries needed, following should work:
```
mkdir openmx-hks
git clone https://github.com/pulkin/openmx-hks.git
cd openmx-hks/src
make
```

## Running examples

Several helpful examples provided:
```
cd openmx-hks/examples/01_hks2json
./run
```

## Bugs

Report [here](https://github.com/pulkin/openmx-hks/issues).
