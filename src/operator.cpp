#include <stdexcept>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/Util/SelectionRule.h>

#include "operator.hpp"

using namespace Spectra;

Eigen::VectorXcd IRLM_eigen(Operator &op, int nb_eigen, Eigen::MatrixXcd& eigenvectors) {
    SparseGenMatProd<double> mat_prod(op);
    GenEigsSolver<SparseGenMatProd<double>> eigs(mat_prod, nb_eigen, 2 * nb_eigen + 1);

    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::SmallestReal);
    if (eigs.info() != Spectra::CompInfo::Successful) {
        throw std::runtime_error("Eigenvalue computation failed.");
    }

    Eigen::VectorXcd eigenvalues = eigs.eigenvalues();
    eigenvectors = eigs.eigenvectors();

    return eigenvalues;
}

Eigen::VectorXd exact_eigen(const Operator &op, Eigen::MatrixXd& eigenvectors) {
    Eigen::MatrixXd dense_smat = Eigen::MatrixXd(op);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(dense_smat);

    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue computation failed.");
    }

    eigenvectors = eigensolver.eigenvectors();
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();

    return eigenvalues;
}

void FOLM_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V) {
    throw std::logic_error("This function has not been implemented yet.");
}

Eigen::VectorXd FOLM_eigen(Operator &op, int nb_iter, Eigen::MatrixXd& eigenvectors) {
    throw std::logic_error("This function has not been implemented yet.");
}

// ============= PHASE TRANSITION CHARACTERIZATION ============= //

// TODO: Implement the order parameter calculation
double order_parameter(const Eigen::VectorXd& eigenvalues, const Eigen::MatrixXd& eigenvectors)  {
    throw std::logic_error("This function has not been implemented yet.");
}

double gap_ratio(Operator &op) {
    Eigen::MatrixXcd eigenvectors;

    int D = op.rows();
    int nb_eigen = std::min(20, D / 2);
    double E0 = std::real(IRLM_eigen(op, nb_eigen, eigenvectors)[0]);
    double E1 = std::real(IRLM_eigen(op, nb_eigen, eigenvectors)[1]);
    return (E1 - E0) / (E1 + E0);
}


/* SPECIFIC CALCULATIONS */


//TODO : A ECRIRE
// void add_chemical_potential(double mu, int n) {
//     *this = *this + Identity(D) * (mu * n);
// }

//TODO : A ECRIRE
// void add_interaction(double U, const Eigen::VectorXd& interaction_matrix) {
//     throw std::logic_error("This function has not been implemented yet.");
// }


///// THERMODYNAMICAL FUNCTIONS /////


double partition_function(const Eigen::VectorXd& eigenvalues, double temperature) {
    const double k_B = 1.380649e-23;
    double beta = 1 / (k_B * temperature);
    return (-beta * eigenvalues).array().exp().sum();
}

// TODO: Implement the canonical density matrix calculation
// void canonical_density_matrix(const Eigen::VectorXd& eigenvalues, double temperature) const {
//     throw std::logic_error("This function has not been implemented yet.");
// }

//TODO : A OPTIMISER
// double boson_density(double dmu, int n) {
//     [[maybe_unused]] Eigen::MatrixXcd eigenvectors;
//     int nb_eigen = std::min(20, D/2);
//     std::complex<double> E0 = this->IRLM_eigen(nb_eigen, eigenvectors)[0]; // ground state energy at mu
//     *this = *this + Identity(D) * dmu* n;
//     std::complex<double> E1 = this->IRLM_eigen(nb_eigen, eigenvectors)[0]; // ground state energy at mu + dmu
//     *this = *this + Identity(D) * (-dmu * n);
//     return -(std::real(E1) - std::real(E0)) / dmu;
// }

//TODO : A OPTIMISER
// double compressibility(double dmu, int n) {
//     double n1 = this->boson_density(dmu, n); // density at mu
//     *this = *this + Identity(D) * dmu * n;
//     double n2 = this->boson_density(dmu, n); // density at mu + dmu
//     *this = *this + Identity(D) * (-dmu * n);
//     return std::pow(n1,-2) * (n2-n1)/dmu;
// }