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

- `docs/`: Contains the generated documentation for the project.

- `CMakeLists.txt`: CMake configuration file for building the project.

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html): A C++ template library for linear algebra.
- [Spectra](https://spectralib.org/): A header-only C++ library for large scale eigenvalue problems.
- [Doxygen](http://www.doxygen.nl/): A tool for generating documentation from annotated C++ sources.
  
## Building the Project

To build the project, follow these steps:

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/Bose-Hubbard-Phase-Transition.git
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

1. Ensure Doxygen is installed on your system.

2. Run Doxygen with the provided `Doxyfile`:
    ```sh
    doxygen Doxyfile
    ```

The documentation will be generated in the `docs` directory.

## Usage

The project provides classes to represent the Bose-Hubbard Hamiltonian and to perform various operations and calculations related to the model. The main classes are:

- `BH`: Represents the Bose-Hubbard Hamiltonian.
- `Neighbours`: Generates the list of neighbours for different lattice structures i.e. the geometry of the lattice.
- `Operator`: Provides various matrix operations and diagonalization methods.

## Authors

- [Maximilien HANTONNE](https://github.com/Maximilien-Hantonne)
- [Remy LYSCAR](https://github.com/Remy-Lyscar)
  
## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
