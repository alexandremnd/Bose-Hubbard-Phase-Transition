#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <omp.h>

#include "hamiltonian.h"
#include "operator.hpp"
#include "neighbours.h"
#include "clock.hpp"


void compute_phase_transition(std::string fixed_param, double fixed, double p1_min, double p1_max, double p2_min, double p2_max, double s, Operator& JH, Operator& UH, Operator& uH, std::ofstream& file) {
    double J = 0;
    double *parameter_1, *parameter_2;

    if (fixed_param == "J") {
        J = fixed;
    }

    size_t size_x = static_cast<int>((p1_max - p1_min) / s) + 1;
    size_t size_y = static_cast<int>((p2_max - p2_min) / s) + 1;
    size_t size = size_x * size_y;
    std::vector<double> x_values(size);
    std::vector<double> y_values(size);
    std::vector<double> gap_values(size);
    std::vector<double> density_values(size);
    std::vector<double> compressibility_values(size);

    Operator H = Eigen::SparseMatrix<double>(JH.size(), JH.size());


    for (size_t i = 0; i < size_x; ++i) {
        for (size_t j = 0; j < size_y; ++j) {
            double U = p1_min + i * s;
            double mu = p2_min + j * s;

            H = JH * J + UH * U + uH * mu;

            double gap_value = gap_ratio(H);
            double boson_value = 0;
            double compress_value = 0;

            x_values[i*size_y + j] = U;
            y_values[i*size_y + j] = mu;
            gap_values[i*size_y + j] = gap_value;
            density_values[i*size_y + j] = boson_value;
            compressibility_values[i*size_y + j] = compress_value;
        }
    }

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            file << x_values[i*size_y + j] << " " << y_values[i*size_y + j] << " " << gap_values[i*size_y + j] << " " << density_values[i*size_y + j] << " " << compressibility_values[i*size_y +j] << "\n";
        }
    }

    return;
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
    Clock c;
    c.start();
    BH jmatrix(nei, m, n, 1, 0, 0);
    Operator JH = jmatrix.getHamiltonian();

    BH Umatrix(nei, m, n, 0, 1, 0);
    Operator UH = Umatrix.getHamiltonian();

    BH umatrix(nei, m, n, 0, 0, 1);
    Operator uH = umatrix.getHamiltonian();


    // CALCULATING AND SAVING THE PHASE TRANSITION PARAMETERS
    if (fixed_param == "J") {
        compute_phase_transition(fixed_param, J, U_min, U_max, mu_min, mu_max, s, JH, UH, uH, file);
    } else if (fixed_param == "U") {
        compute_phase_transition(fixed_param, U, J_min, J_max, mu_min, mu_max, s, JH, UH, uH, file);
    } else if (fixed_param == "u") {
        compute_phase_transition(fixed_param, mu, J_min, J_max, U_min, U_max, s, JH, UH, uH, file);
    }
    file.close();

    c.time_s("Phase transition calculation:");

    // PLOT OF THE PHASE TRANSITION
	int result = system("python3 plot.py");
	if (result != 0) {
		std::cerr << "Error when executing Python script." << std::endl;
		return 1;
	}

	return 0;
}
