# Bose-Hubbard Phase Transition

This repository contains the implementation of the Bose-Hubbard model, which describes interacting bosons on a lattice. The aim of this project is to study the phase transition in the Bose-Hubbard model using various numerical methods.

## Project Structure

- `include/`: Contains the header files for the project.
  - `hamiltonian.h`: Defines the `BH` class representing the Bose-Hubbard Hamiltonian.
  - `neighbours.h`: Defines the `Neighbours` class for generating the list of neighbours for different lattice structures.
  - `operator.h`: Defines the `Operator` class for various matrix operations and diagonalization methods.

- `src/`: Contains the source files for the project.

- `external/`: Contains external dependencies.
  - `spectra/`: Header-only C++ library for large scale eigenvalue problems.

- `documentation.pdf`: Contains the generated documentation for the project.

- `CMakeLists.txt`: CMake configuration file for building the project.

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html): A C++ template library for linear algebra.
- [Spectra](https://spectralib.org/): A header-only C++ library for large scale eigenvalue problems.
- [Doxygen](http://www.doxygen.nl/): A tool for generating documentation from annotated C++ sources.
  
## Building the Project

To build the project, follow these steps:

1. Clone the repository:
    ```sh
    git clone https://github.com/Maximilien-Hantonne/Bose-Hubbard-Phase-Transition.git
    cd Bose-Hubbard-Phase-Transition
    ```

2. Create a build directory and navigate to it:
    ```sh
    mkdir build
    cd build
    ```

3. Run CMake to configure the project:
    ```sh
    cmake ..
    ```

4. Build the project:
    ```sh
    make
    ```

## Generating Documentation

To generate the documentation using Doxygen, follow these steps:

### HTML Documentation

To generate the HTML documentation using Doxygen, follow these steps:

1. Ensure Doxygen is installed on your system.

2. Run Doxygen with the provided `Doxyfile`:
    ```sh
    doxygen Doxyfile
    ```

The documentation will be generated in the `docs/html` directory.

### LaTeX Documentation

To generate the LaTeX documentation using Doxygen, follow these steps:

1. Ensure Doxygen and LaTeX are installed on your system. You can install LaTeX using:
    ```sh
    sudo apt-get install texlive-full
    ```

2. Run Doxygen with the provided `Doxyfile`:
    ```sh
    doxygen Doxyfile
    ```

3. Navigate to the `docs/latex` directory:
    ```sh
    cd docs/latex
    ```

4. Compile the LaTeX documentation:
    ```sh
    make
    ```

The PDF documentation will be generated in the `docs/latex` directory.

## Usage

The project provides classes to represent the Bose-Hubbard Hamiltonian and to perform various operations and calculations related to the model. The main classes are:

- `BH`: Represents the Bose-Hubbard Hamiltonian.
- `Neighbours`: Generates the list of neighbours for different lattice structures i.e. the geometry of the lattice.
- `Operator`: Provides various matrix operations and diagonalization methods.

### Command-Line Options

The `main` function parses command-line arguments to set up the parameters for the Bose-Hubbard model.

Command-line options:
- `-m, --sites`: Number of sites in the lattice.
- `-n, --bosons`: Number of bosons in the lattice.
- `-J, --hopping`: Hopping parameter.
- `-U, --interaction`: On-site interaction parameter.
- `-u, --potential`: Chemical potential.
- `-r, --range`: Range for chemical potential and interaction.
- `-s, --step`: Step for chemical potential and interaction.
- `-h, --help`: Display usage information.


## Authors

- [Maximilien HANTONNE](https://github.com/Maximilien-Hantonne)
- [Rémy LYSCAR](https://github.com/Remy-Lyscar)
  
## Contributing

Contributions are currently not open :(

## License

This project is licensed under the MIT License. See the LICENSE file for details.
