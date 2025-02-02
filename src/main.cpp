#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>
#include <complex>
#include <getopt.h>

#include "hamiltonian.h"
#include "operator.h"
#include "neighbours.h"


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

    /// PARAMETERS OF THE MODEL
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

	/// GEOMETRY OF THE LATTICE
	Neighbours neighbours(m);
	neighbours.chain_neighbours(); // list of neighbours
	const std::vector<std::vector<int>>& nei = neighbours.getNeighbours();

    // PLOT OF THE PHASE TRANSITION
    double dmu = 0.05 * std::min({J > 0 ? J : std::numeric_limits<double>::max(),
                                  U > 0 ? U : std::numeric_limits<double>::max(),
                                  mu > 0 ? mu : std::numeric_limits<double>::max()});

    std::ofstream file("phase.txt");
    file << fixed_param << " "; 
    if (fixed_param == "J") {
        file << J << std::endl;
    } else if (fixed_param == "U") {
        file << U << std::endl;
    } else if (fixed_param == "u") {
        file << mu << std::endl;
    } else {
        std::cerr << "Error: Invalid fixed parameter specified.\n";
        print_usage();
        return 1;
    } 
    
    if (fixed_param == "J") {
        for (double U = U_min; U <= U_max; U += s) {
            for (double mu = mu_min; mu <= mu_max; mu += s) {
                BH hmatrix(nei, m, n, J, U, mu);
                Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
                Operator H(std::move(smatrix));
                double gap_ratio = H.gap_ratio();
                double boson_density = H.boson_density(dmu, n);
                double compressibility = H.compressibility(dmu, n);
                file << U << " " << mu << " " << gap_ratio << " " << boson_density << " " << compressibility << std::endl;
            }
        }
    } else if (fixed_param == "U") {
        for (double J = J_min; J <= J_max; J += s) {
            for (double mu = mu_min; mu <= mu_max; mu += s) {
                BH hmatrix(nei, m, n, J, U, mu);
                Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
                Operator H(std::move(smatrix));
                double gap_ratio = H.gap_ratio();
                double boson_density = H.boson_density(dmu, n);
                double compressibility = H.compressibility(dmu, n);
                file << J << " " << mu << " " << gap_ratio << " " << boson_density << " " << compressibility << std::endl;
            }
        }
    } else if (fixed_param == "u") {
        for (double J = J_min; J <= J_max; J += s) {
            for (double U = U_min; U <= U_max; U += s) {
                BH hmatrix(nei, m, n, J, U, mu);
                Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
                Operator H(std::move(smatrix));
                double gap_ratio = H.gap_ratio();
                double boson_density = H.boson_density(dmu, n);
                double compressibility = H.compressibility(dmu, n);
                file << J << " " << U << " " << gap_ratio << " " << boson_density << " " << compressibility << std::endl;
            }
        }
    } 
    file.close();
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
