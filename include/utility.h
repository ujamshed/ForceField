#pragma once
#include <armadillo>
#include <cmath>

/**
 * @brief Function to calculate the difference between 2 vectors in 3D space.
 * 
 * @param atom_i (arma::rowvec): x, y, z coordinate information for atom i.
 * @param atom_j (arma::rowvec): x, y, z coordinate information for atom j.
 * @return arma::rowvec 
 */
arma::rowvec vector_subtraction(arma::rowvec atom_i, arma::rowvec atom_j);

/**
 * @brief Function to calculate the euclidian distance between 2 atoms in 3D space.
 * 
 * @param atom_i (arma::rowvec): x, y, z coordinate information for atom i.
 * @param atom_j (arma::rowvec): x, y, z coordinate information for atom j.
 * @return double 
 */
double euclidian_distance(arma::rowvec atom_i, arma::rowvec atom_j);

/**
 * @brief Function to calculate the angle between 3 atoms in 3D space (in degrees). Using law of cosines.
 * 
 * @param atom_i (arma::rowvec): x, y, z coordinate information for atom i.
 * @param atom_j (arma::rowvec): x, y, z coordinate information for atom j.
 * @param atom_k (arma::rowvec): x, y, z coordinate information for atom k.
 * @return double 
 */
double calc_angle(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k);

/**
 * @brief Function to calculate the torsional angle between 4 atoms in 3D space (in degrees).
 * 
 * @param atom_i (arma::rowvec): x, y, z coordinate information for atom i.
 * @param atom_j (arma::rowvec): x, y, z coordinate information for atom j.
 * @param atom_k (arma::rowvec): x, y, z coordinate information for atom k.
 * @param atom_l (arma::rowvec): x, y, z coordinate information for atom l.
 * @return double 
 */
double calc_dihedral(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k, arma::rowvec atom_l);

/**
 * NOT IMPLEMENTED 
 * 
 * Helper function to calculate the partial atomic charge between an atom and its bonded atoms.
 * Need to send in the information about an atom and ONLY its bonds.
 */
double calc_partial_atomic_charge(arma::mat atom_bonds);