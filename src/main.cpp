#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include <complex>
#include <getopt.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <omp.h>

#include "hamiltonian.h"
#include "operator.h"
#include "neighbours.h"

void compute_gap_ratio(std::string fixed_param, double J, double U, double mu, Operator& JH, Operator& UH, Operator& uH, std::ofstream& file) {
    Operator H = JH * J + UH * U + uH * mu;

    double gap_ratio = H.gap_ratio();
    double boson_density = 0;
    double compressibility = 0;

    double p1, p2;
    if (fixed_param == "J") {
        p1 = U;
        p2 = mu;
    } else if (fixed_param == "U") {
        p1 = J;
        p2 = mu;
    } else if (fixed_param == "u") {
        p1 = J;
        p2 = U;
    } else {
        std::cerr << "Error: Invalid fixed parameter specified.\n";
        return;
    }

    std::stringstream ss;
    ss << p1 << " " << p2 << " " << gap_ratio << " " << boson_density << " " << compressibility << "\n";
    file << ss.str();
}

void compute_phase_transition(std::string fixed_param, double fixed, double p1_min, double p1_max, double p2_min, double p2_max, double s, Operator& JH, Operator& UH, Operator& uH, std::ofstream& file) {
    double J = 0, U = 0, mu = 0;
    double *parameter_1, *parameter_2;

    if (fixed_param == "J") {
        J = fixed;
        parameter_1 = &U;
        parameter_2 = &mu;
    } else if (fixed_param == "U") {
        U = fixed;
        parameter_1 = &J;
        parameter_2 = &mu;
    } else if (fixed_param == "u") {
        mu = fixed;
        parameter_1 = &J;
        parameter_2 = &U;
    } else {
        std::cerr << "Error: Invalid fixed parameter specified.\n";
        return;
    }

    for (int i = 0; i < static_cast<int>((p1_max - p1_min) / s) + 1; ++i) {
        for (int j = 0; j < static_cast<int>((p2_max - p2_min) / s) + 1; ++j) {
            *parameter_1 = p1_min + i * s;
            *parameter_2 = p2_min + j * s;

            compute_gap_ratio(fixed_param, J, U, mu, JH, UH, uH, file);
        }
    }

    return;
}

/**
 * @brief Get the memory usage of the program in KB and optionally print it.
 * @param print If true, print the memory usage.
 * @return The memory usage in KB.
 */
long get_memory_usage(bool print = false) {
    std::ifstream statm_file("/proc/self/statm");
    long memory_usage = -1;

    if (statm_file.is_open()) {
        long size, resident, share, text, lib, data, dt;
        statm_file >> size >> resident >> share >> text >> lib >> data >> dt;
        memory_usage = resident * (sysconf(_SC_PAGESIZE) / 1024); // Convert pages to KB
    }

    if (memory_usage == -1) {
        std::cerr << "Error reading memory usage from /proc/self/statm." << std::endl;
    } else if (print) {
        std::cout << "Memory usage: " << memory_usage << " KB" << std::endl;
    }

    return memory_usage;
}


/**
 * @brief Get the available memory in KB.
 */
long get_available_memory() {
    struct sysinfo info;
    if (sysinfo(&info) != 0) {
        std::cerr << "Error getting system info." << std::endl;
        return -1;
    }
    return info.freeram / 1024;
}


size_t estimateSparseMatrixMemoryUsage(const Eigen::SparseMatrix<double>& matrix) {
    size_t numNonZeros = matrix.nonZeros();
    size_t numCols = matrix.cols();
    size_t memoryUsage = numNonZeros * sizeof(double);
    memoryUsage += numNonZeros * sizeof(int);
    memoryUsage += (numCols + 1) * sizeof(int);
    return memoryUsage;
}

/**
 * @brief Print the usage information for the program.
 */
void print_usage() {
    std::cout << "Usage: program [options]\n"
              << "Options:\n"
              << "  -m, --sites       Number of sites\n"
              << "  -n, --bosons      Number of bosons\n"
              << "  -J, --hopping     Hopping parameter\n"
              << "  -U, --interaction On-site interaction\n"
              << "  -u, --potential   Chemical potential\n"
              << "  -r, --range     Range for varying parameters\n"
              << "  -s, --step      Step for varying parameters (with s < r)\n"
              << "  -f --fixed      Fixed parameter (J, U or u) \n";
}


/**
 * @brief Main function for the Bose-Hubbard Phase Transition program.
 *
 * This function parses command-line arguments to set up the parameters for the Bose-Hubbard model.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 *
 * Command-line options:
 * - `-m, --sites`: Number of sites in the lattice.
 * - `-n, --bosons`: Number of bosons in the lattice.
 * - `-J, --hopping`: Hopping parameter.
 * - `-U, --interaction`: On-site interaction parameter.
 * - `-u, --potential`: Chemical potential.
 * - `-r, --range`: Range for chemical potential and interaction.
 * - `-s, --step`: Step for chemical potential and interaction.
 * - `-f, --fixed`: Fixed parameter (J, U, or u).
 * - `-h, --help`: Display usage information.
 *
 * @return int Exit status of the program.
 * @warning s must be smaller than r.
 */

int main(int argc, char *argv[]) {

    // // INITIAL MEMORY USAGE
    // long base_memory = get_memory_usage();

    // PARAMETERS OF THE MODEL
    int m, n;
    double J, U, mu, s, r;
    [[maybe_unused]] double J_min, J_max, mu_min, mu_max, U_min, U_max;
    std::string fixed_param;

    const char* const short_opts = "m:n:J:U:u:r:s:f:h";
    const option long_opts[] = {
        {"sites", required_argument, nullptr, 'm'},
        {"bosons", required_argument, nullptr, 'n'},
        {"hopping", required_argument, nullptr, 'J'},
        {"interaction", required_argument, nullptr, 'U'},
        {"potential", required_argument, nullptr, 'u'},
        {"range", required_argument, nullptr, 'r'},
        {"step", required_argument, nullptr, 's'},
        {"fixed", required_argument, nullptr, 'f'},
		{"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
        if (-1 == opt) break;
        switch (opt) {
            case 'm':
                m = std::stoi(optarg);
                break;
            case 'n':
                n = std::stoi(optarg);
                break;
            case 'J':
                J = std::stod(optarg);
                break;
            case 'U':
                U = std::stod(optarg);
                break;
            case 'u':
                mu = std::stod(optarg);
                break;
            case 'r':
                r = std::stod(optarg);
                break;
            case 's':
                s = std::stod(optarg);
                break;
            case 'f':
                fixed_param = optarg;
                break;
			case 'h':
            default:
                print_usage();
                return 0;
        }
    }

    J_min = J;
    J_max = J + r;
	mu_min = mu;
    mu_max = mu + r;
    U_min = U;
    U_max = U + r;

    [[maybe_unused]] double dmu = 0.05 * std::min({J > 0 ? J : std::numeric_limits<double>::max(),
                                  U > 0 ? U : std::numeric_limits<double>::max(),
                                  mu > 0 ? mu : std::numeric_limits<double>::max(),
                                  s > 0 ? s : std::numeric_limits<double>::max()});


	// GEOMETRY OF THE LATTICE
	Neighbours neighbours(m);
	neighbours.chain_neighbours(); // list of neighbours
	const std::vector<std::vector<int>>& nei = neighbours.getNeighbours();


    // OPENING THE FILE TO SAVE THE MAP VALUES
    std::ofstream file("phase.txt");
    file << fixed_param << " ";
    if (fixed_param == "J") {
        if (J == 0) {
            std::cerr << "Error: Fixed parameter J cannot be zero.\n";
            return 1;
        }
        file << J << std::endl;
    } else if (fixed_param == "U") {
        if (U == 0) {
            std::cerr << "Error: Fixed parameter U cannot be zero.\n";
            return 1;
        }
        file << U << std::endl;
    } else if (fixed_param == "u") {
        if (mu == 0) {
            std::cerr << "Error: Fixed parameter mu cannot be zero.\n";
            return 1;
        }
        file << mu << std::endl;
    } else {
        std::cerr << "Error: Invalid fixed parameter specified.\n";
        print_usage();
        return 1;
    }


    // HAMILTONIAN INITIALIZATION
    auto start = std::chrono::high_resolution_clock::now();
    BH jmatrix(nei, m, n, 1, 0, 0);
    Eigen::SparseMatrix<double> jsmatrix = jmatrix.getHamiltonian();
    Operator JH(std::move(jsmatrix));

    BH Umatrix(nei, m, n, 0, 1, 0);
    Eigen::SparseMatrix<double> Usmatrix = Umatrix.getHamiltonian();
    Operator UH(std::move(Usmatrix));

    BH umatrix(nei, m, n, 0, 0, 1);
    Eigen::SparseMatrix<double> usmatrix = umatrix.getHamiltonian();
    Operator uH(std::move(usmatrix));


    // SETTING THE NUMBER OF THREADS FOR PARALLELIZATION
    long available_memory = get_available_memory();
    std::cout << "Available memory: " << available_memory << " KB" << std::endl;
    size_t jsmatrixMemoryUsage = estimateSparseMatrixMemoryUsage(jsmatrix);
    size_t usmatrixMemoryUsage = estimateSparseMatrixMemoryUsage(Usmatrix);
    size_t umatrixMemoryUsage = estimateSparseMatrixMemoryUsage(usmatrix);
    size_t totalMemoryUsage = jsmatrixMemoryUsage + usmatrixMemoryUsage + umatrixMemoryUsage;
    std::cout << "Estimated total memory usage for sparse matrices: " << totalMemoryUsage << " bytes" << std::endl;
    int num_threads = std::min(available_memory / totalMemoryUsage, static_cast<size_t>(omp_get_max_threads()));
    omp_set_num_threads(8);
    std::cout << "Using OpenMP with " << num_threads << " threads." << std::endl;


    // CALCULATING AND SAVING THE PHASE TRANSITION PARAMETERS
    if (fixed_param == "J") {
        compute_phase_transition(fixed_param, J, U_min, U_max, mu_min, mu_max, s, JH, UH, uH, file);
    } else if (fixed_param == "U") {
        compute_phase_transition(fixed_param, U, J_min, J_max, mu_min, mu_max, s, JH, UH, uH, file);
    } else if (fixed_param == "u") {
        compute_phase_transition(fixed_param, mu, J_min, J_max, U_min, U_max, s, JH, UH, uH, file);
    }
    file.close();


    // EFFICIENCY OF THE CALCULATION
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    if (duration.count() > 60) {
        int minutes = static_cast<int>(duration.count()) / 60;
        double seconds = duration.count() - (minutes * 60);
        std::cout << "Parameters calculation time: " << minutes << " minutes " << seconds << " seconds" << std::endl;
    } else {
        std::cout << "Parameters calculation time: " << duration.count() << " seconds" << std::endl;
    }
    get_memory_usage(true);


    // PLOT OF THE PHASE TRANSITION
	int result = system("python3 plot.py");
	if (result != 0) {
		std::cerr << "Error when executing Python script." << std::endl;
		return 1;
	}

	// /// DIAGONALIZATION OF THE HAMILTONIAN
	// BH hmatrix (nei, m, n, J, U, mu);
	// Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
	// Operator H(std::move(smatrix));

	// // USING THE FOLM
	// int k = H.size(); // Number of eigenvalues to calculate
	// Eigen::VectorXd v_0 = Eigen::VectorXd::Random(H.size()); // Random initial vector
	// Eigen::MatrixXd eigenvectors1;
	// Eigen::MatrixXd V;
	// auto start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXd eigenvalues1 = H.FOLM_eigen(k, v_0,eigenvectors1); // FOLM
	// auto end = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> duration = end - start;
	// std::cout << "FOLM execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << eigenvalues1.transpose()[0] << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues1.size() << std::endl << std::endl;

	// // USING THE IRLM
    // Eigen::MatrixXd eigenvectors2;
	// start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXcd eigenvalues2 = H.IRLM_eigen(1, eigenvectors2); // IRLM
	// end = std::chrono::high_resolution_clock::now();
	// duration = end - start;
	// std::cout << "IRLM execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << std::real(eigenvalues2.transpose()[0]) << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues2.size() << std::endl << std::endl;

	// // USING EXACT DIAGONALIZATION
	// Eigen::MatrixXd eigenvectors3;
	// start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXd eigenvalues3 = H.exact_eigen(eigenvectors3); // Exact diagonalization
	// end = std::chrono::high_resolution_clock::now();
	// duration = end - start;
	// std::cout << "Exact diagonalization execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << eigenvalues3.transpose()[0] << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues3.size() << std::endl << std::endl;


	/// PHASE TRANSITION CALCULATIONS
	// double boson_density = H.boson_density(0.1, n);
	// std::cout << "boson density : " << boson_density << std::endl;

	// double compressibility = H.compressibility(0.1, n);
	// std::cout << "isothermal compressibility : " << compressibility << std::endl << std::endl;

	return 0;
}
