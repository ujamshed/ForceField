#pragma once
#include "input.h"
#include "utility.h"
#include "energy.h"

// Function to calculate central difference using input matrix
arma::mat grad(arma::mat Coordinates, double step_size, int num_atoms, calcEnergy& Molecule);

// Steepest Descent
// Trivial, will not work for all starting arrangements
arma::mat steepest_descent(arma::mat Coordinates, double step_size, double tol, int num_atoms, calcEnergy& Molecule);

// Backtrack line-search with wolfe conditions
double line_search(int num_atoms, arma::mat coordinates, arma::mat search_direction, arma::mat gradient, calcEnergy& Molecule);

// BFGS
arma::mat BFGS(arma::mat Coordinates, double tol, int num_atoms, calcEnergy& Molecule);