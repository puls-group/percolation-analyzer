# Percolation Detection Analysis

## Description

In a finite MD system with periodic boundary conditions (pbc) of any (triclinic) shape, it is relevant to detect the percolation point in order to formally define the point of gelation. 
Here the percolation point is defined as the time, at which a molecule first percolates (= connects to itself) in each dimension. 

This repository contains the necessary code to analyze whether or not a molecule is percolating, in how many dimensions it does so (1d -> string like, 2d -> sheet like, 3d -> grid-like) and returns the information for each individual connected molecule in the graph provided.

## License/Permission of use

The code is provided as-is without any warranty for its correctness under the MIT license (see `LICENSE` file). We have tested it and are reasonably sure that it works, but it is always possible for bugs to occur in special cases. If you think you have found a bug, please let us know via the issues of this project.

If you use the code provided in this repository for research purposes, please cite the paper "An Exact Algorithm to Detect the Percolation Transition in Molecular Dynamics Simulations of Cross-Linking Polymer Networks" (doi: [10.1021/acs.jctc.1c00423 ](https://doi.org/10.1021/acs.jctc.1c00423 )) as a reference.

BibTeX:
```
@article{Livraghi_2021, doi = {10.1021/acs.jctc.1c00423}, url = {https://doi.org/10.1021%2Facs.jctc.1c00423}, year = 2021, month = {sep}, publisher = {American Chemical Society ({ACS})}, author = {Mattia Livraghi and Kevin H\"ollring and Christian R. Wick and David M. Smith and Ana-Sun{\v{c}}ana Smith}, title = {An Exact Algorithm to Detect the Percolation Transition in Molecular Dynamics Simulations of Cross-Linking Polymer Networks}, journal = {Journal of Chemical Theory and Computation} } 
```

This repository also provides a doi for each full release of the software. Please consult the releases page for the full list of releases.

## How to use

### General Percolation Graph interface
We will explain the workings of this library by example of the cpp classes provided in `include/percolation-detection.hpp` and `include/molecular-graph.hpp`. There is also a C-wrapper for the percolation backend provided in `include/percolation-analyzer.h`.

The design decisions for the interface are the same.

In order to use the percolation analysis, you need to build a `percolation::PercolationGraph` object for each frame of a trajectory and run the method `percolation::PercolationGraph::get_component_percolation_info()` method on it, which returns the percolation information for each of the molecules within the frame as a `std::vector`. 
The `percolation_dim` member variable of the vector entries contains the percolation dimension (0: no percolation, 1: line structure, 2: sheet structure, 3: full grid structure) while the `vertices` member holds the data (including index) associated with the atoms/vertices within the respective molecule (for identification purposes).

To build the `percolation::PercolationGraph` object, you need to declare the number of atoms in your system in `percolation::PercolationGraph::reserve_vertices()` to reserve the memory and then register the vertex data in `percolation::PercolationGraph::add_vertex()` and the link/bond/edge data in `percolation::PercolationGraph::add_edge()`.
If you want to include more information than currently provided by the library, you can extend the structures `percolation::VertexData` and `percolation::EdgeData`  to allow for more data being passed in and out of the analysis.

Each edge needs to be provided with an integer translation vector as a `percolation::TranslationVector` object, to detail the pbc crossings as detailed in our publication explaining the percolation detection algorithm (entry +1 if pbc crossed upwards from source to head, -1 of pbc crossed downwards from source to head, 0 if edge completely within pbc cell).

### The Molecular Graph interface

To simplify the building of the `percolation::PercolationGraph` object, we provide a helper class `mol::MolecularGraph` in `include/molecular-graph.hpp` in which you can simply provide the pbc information as a triclinic base via `mol::MolecularGraph::set_basis()`, the information for each atom/vertex via `mol::MolecularGraph::set_atom_position()` and the bond information via `mol::MolecularGraph::add_bond()` which only takes the information, which atoms are bonded. Please take note, that the MolecularGraph class only converts to the PercolationGraph class correctly, if all bonds only ever cross over in up to one of the next neighboring pbc cells. If your bonds may cross one full pbc cell or more, you need to build the PercolationGraph yourself.

To convert the `mol::MolecularGraph` object into the corresponding `percolation::PercolationGraph` object, you need to call `mol::MolecularGraph::get_percolation_graph()`. The source code of that function in `src/molecular-graph.cpp` can also be used as an illustration on how to convert the molecular graph into the percolation graph in general.

### Building the library/Build system

The library provides a build system based on cmake (so you will need to install that before attempting a build of the repository). 
In order to build all components, the library "Catch2" is used for testing, so if you want to build all targets, please include all submodules after cloning using the command `git submodule update --init` after your `git clone` command to download the Catch2-files.

To build the library and sample program provided in `src/sample_graph_builder.cpp` we suggest running the following commands if you are not familiar with the cmake build system:
```
mkdir build
cd build 
cmake ..
```

This will set up the cmake configuration in a subdirectory `build` of the current directory.
* To build the cpp percolation library (output: `bin/percolation-analyzer-cpp.lib`) you can then run `make percolation-analyzer-cpp` in the `build` directory
* To build the c percolation library wrapper (output: `bin/percolation-analyzer-c.lib`) you can then run `make percolation-analyzer-c` in the `build` directory
* To build the cpp percolation library with the molecular graph helper class (output: `bin/molecular-graph-cpp.lib`) you can then run `make molecular-graph-cpp` in the `build` directory
* To build the sample program setting up a random graph, (output: `bin/sample_graph`) you can then run `make sample_graph` in the `build` directory
* To build the tests to check for the correct operation of the percolation analysis (output: `bin/test_runner`) you can then run `make test_runner` in the `build` directory

In order for your own program to use this library, you need to link against the appropriate library that you want to use ( `bin/percolation-analyzer-cpp.lib`,  `bin/percolation-analyzer-c.lib` or  `bin/molecular-graph-cpp.lib`) as well as add the files in the `include` directory to your include path.
## Authors

The code was written by Kevin Höllring in 2021, phd student at PULS group of Friedrich-Alexander-University Erlangen-Nuremberg.
