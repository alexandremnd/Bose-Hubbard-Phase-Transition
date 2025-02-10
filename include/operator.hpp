#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

typedef Eigen::SparseMatrix<double> Operator;


// ==========================
//     DIAGONALIZATION
// ==========================

/**
 * @brief Implicitly Restarted Lanczos Method (IRLM) for a sparse matrix
 *
 * @param op Operator to diagonalize
 * @param nb_eigen Number of eigenvalues to compute
 * @param eigenvectors Returned eigenvectors
 * @return Eigen::VectorXcd Eigenvalues
 */
Eigen::VectorXcd IRLM_eigen(Operator &op, int nb_eigen, Eigen::MatrixXcd& eigenvectors);

Eigen::VectorXd IRLM_eigen_sym(Operator &op, int nb_eigen, Eigen::MatrixXd& eigenvectors);

Eigen::VectorXd IRLM_eigen_herm(Operator &op, int nb_eigen, Eigen::MatrixXd& eigenvectors);

/**
* @brief Calculate the exact eigenvalues and eigenvectors of the Hamiltonian by an exact diagonalization.
*
* @param op Operator to diagonalize
* @param eigenvectors An empty matrix to store the eigenvectors.
* @return Eigen::Matrix<double> The vector of eigenvalues.
* @warning This function is computationally expensive and should be used with caution.
*/
Eigen::VectorXd exact_eigen(const Operator &op, Eigen::MatrixXd& eigenvectors);

/**
 * @brief Full Orthogonalization Lanczos Method (FOLM) for a sparse matrix
 *
 * @param op Operator to diagonalize
 * @param nb_iter Number of iterations of the algorithm
 * @param v_0 First vector of the Krylov basis
 * @param T
 * @param V
 */
void FOLM_diag(Operator &op, int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V);

/**
* @brief Calculate the approximate eigenvalues and eigenvectors of the Hamiltonian using the Full Orthogonalization Lanczos Method.
*
* @param op Operator to diagonalize
* @param nb_iter The number of iterations.
* @param eigenvectors An empty matrix to store the eigenvectors.
* @return Eigen::VectorXd The vector of eigenvalues.
* @warning Ensure that nb_iter is greater than 1. The calculation might be wrong if the number of iterations is too low.
*/
Eigen::VectorXd FOLM_eigen(Operator &op, int nb_iter, Eigen::MatrixXd& eigenvectors);

/**
 * @brief Sorts the eigenvalues and eigenvectors in descending order
 *
 * @param eigenvalues Eigenvalues to order
 * @param eigenvectors Eigenvectors to order
 */
void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors);


// PHASE TRANSITION CALCULATIONS :

/**
* @brief Calculate the order parameter of the system.
*
* @param eigenvalues The vector of eigenvalues.
* @param eigenvectors The matrix of eigenvectors.
* @return double The order parameter.
*/
double order_parameter(const Eigen::VectorXd& eigenvalues, const Eigen::MatrixXd& eigenvectors);

/**
* @brief Calculate the energy gap ratio of the system.
*
* @param eigenvalues The vector of eigenvalues.
* @return double The energy gap ratio.
*/
double gap_ratio(Operator &op);

/**
* @brief Add a potential term to the operator.
*/
void add_chemical_potential(double mu, int n);

/**
* @brief Add an interaction term to the operator.
*/
void add_interaction(double U, const Eigen::VectorXd& basis);

// ===========================
// THERMODYNAMICAL FUNCTIONS :

/**
* @brief Calculate the partition function Z for an already diagonalized Hamiltonian.
*
* @param eigenvalues The vector of eigenvalues.
* @param temperature The temperature.
* @return double The partition function.
*/
double partition_function(const Eigen::VectorXd& eigenvalues, double temperature);

/**
* @brief Calculate the canonical density matrix for an already diagonalized Hamiltonian.
*
* @param eigenvalues The vector of eigenvalues.
* @param temperature The temperature.
*/
void canonical_density_matrix(const Eigen::VectorXd& eigenvalues, double temperature);

/**
* @brief Calculate the boson density of the system.
*
* @param dmu The difference of chemical potential.
* @param n The number of bosons.
* @return double The boson density.
* @warning This function modifies the values of the operator from mu to mu + dmu.
*/
double boson_density(double dmu, int n);

/**
* @brief Calculate the isothermal compressibility of the system.
*
* @param dmu The difference of chemical potential.
* @param n The number of bosons.
* @return double The compressibility.
* @warning This function modifies the values of the operator from mu to mu + dmu.
*/
double compressibility(double dmu, int n);