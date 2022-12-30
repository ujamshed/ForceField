#pragma once
#include "input.h"
#include "utility.h"
#include "energy.h"

/**
 * @brief Function to calculate central difference using input matrix
 * 
 * @param Coordinates (arma::mat): Coordinates of the molecule.
 * @param step_size (double): Step size to adjust by when computing the central difference.
 * @param num_atoms (int): Number of atoms in the molecule.
 * @param Molecule (calcEnergy): calcEnergy object which contains the molecule information and can compute the total molecular energy.
 * @return arma::mat 
 */
arma::mat grad(arma::mat Coordinates, double step_size, int num_atoms, calcEnergy& Molecule);

/**
 * @brief Steepest descent algorithm
 * @details Trivial optimization method and will not work for most molecules.
 * 
 * @param Coordinates (arma::mat): Coordinates of the molecule.
 * @param step_size (double): Step size to adjust by when computing the central difference.
 * @param tol (double): Tolerance after which to stop iterating.
 * @param num_atoms (int): Number of atoms in the molecule.
 * @param Molecule (calcEnergy): calcEnergy object which contains the molecule information and can compute the total molecular energy.
 * @return arma::mat 
 */
arma::mat steepest_descent(arma::mat Coordinates, double step_size, double tol, int num_atoms, calcEnergy& Molecule);

/**
 * @brief Backtrack line-search with wolfe conditions
 * 
 * @param num_atoms (int): Number of atoms in the molecule.
 * @param coordinates (arma::mat): Coordinates of the molecule.
 * @param search_direction (arma::mat) Direction of where to start optimizing a.
 * @param gradient (arma::mat) Gradient of the current coordinate system.
 * @param Molecule (calcEnergy): calcEnergy object which contains the molecule information and can compute the total molecular energy.
 * @return double 
 */
double line_search(int num_atoms, arma::mat coordinates, arma::mat search_direction, arma::mat gradient, calcEnergy& Molecule);

/**
 * @brief BFGS algorithm
 * 
 * @param Coordinates (arma::mat): Coordinates of the molecule.
 * @param tol (double): Tolerance after which to stop iterating.
 * @param num_atoms (int): Number of atoms in the molecule.
 * @param Molecule (calcEnergy): calcEnergy object which contains the molecule information and can compute the total molecular energy.
 * @return arma::mat 
 */
arma::mat BFGS(arma::mat Coordinates, double tol, int num_atoms, calcEnergy& Molecule);