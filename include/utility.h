#pragma once
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>

// Helper function to calculate the difference between 2 vectors in 3D space.
arma::rowvec vector_subtraction(arma::rowvec atom_i, arma::rowvec atom_j);

// Helper function to calculate the euclidian distance between 2 atoms in 3D space.
double euclidian_distance(arma::rowvec atom_i, arma::rowvec atom_j);

// Helper function to calculate the angle between 3 atoms in 3D space (in degrees).
// Using law of cosines
double calc_angle(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k);

// Helper function to calculate the torsional angle between 4 atoms in 3D space (in degrees).
double calc_dihedral(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k, arma::rowvec atom_l);

// Helper function to calculate the partial atomic charge between an atom and its bonded atoms.
// Need to send in the information about an atom and ONLY its bonds.
double calc_partial_atomic_charge(arma::mat atom_bonds);

// Helper function to output the optimized coordinates into an sdf file to view in pymol or related molecule viewer.
void sdf_output(int num_atoms, arma::irowvec atom_identity, arma::mat coordinates, arma::imat bonding, std::string name);
