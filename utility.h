#include <armadillo>
#include <cmath>

// Helper function to calculate the euclidian distance between 2 atoms in 3D space.
double euclidian_distance(arma::rowvec atom_i, arma::rowvec atom_j);

// Helper function to calculate the angle between 3 atoms in 3D space (in degrees).
// Using law of cosines
double calc_angle(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k);
