#ifndef OPERATOR_H
#define OPERATOR_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <stdexcept>

class Operator {
private:
    Eigen::SparseMatrix<double> smat;
    double ref; // threshold under which a value is considered null

    /* implement the Lanczos algorithm for a sparse matrix for nb_iter iterations starting with vector v_0 */
    void lanczos_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V);

    void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors); // sort eigenvalues and eigenvectors in descending order

public:
    Operator(int D); // A CHANGER

    /* add a matrix to an operand of type SparseMatrix with same size */
    Operator operator + (const Hamiltonian& operand);

    /* multiply a sparse matrix by a multiplicand of type SparseMatrix with same size */
    Operator operator * (const Hamiltonian& multiplicand);

    /* multiply a sparse matrix by a vector with concomitant size */
    Eigen::VectorXd operator * (const Eigen::VectorXd& vector);

    /* calculate the approximate eigenvalues and eigenvectors of the Hamiltonian in Krylov space using the Lanczos algorithm */
    Eigen::VectorXd lanczos_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V);

    /* calculate the exact eigenvalues and eigenvectors of the Hamiltonian by an exact diagonalization */
    Eigen::VectorXd exact_eigen(Eigen::MatrixXd& eigenvectors);
};

#endif // HAMILTONIAN_H
