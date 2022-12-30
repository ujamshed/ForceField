#pragma once
#include <armadillo>
#include <fstream>
#include <iomanip>

/**
 * @brief Function to output the optimized coordinates into an sdf file to view in pymol or related molecule viewer.
 * 
 * @param num_atoms (int): Number of atoms in the molecule.
 * @param atom_identity (arma::irowvec): Row vector of integers with the identity of the atoms.
 * @param coordinates  (arma::mat): Coordinates of the molecule.
 * @param bonding  (arma::imat): Bonding information of the molecule.
 * @param name (string): Output file name. File extension will be .sdf.
 * @param energy (double): Energy of the current coordinate configuration.
 */
void sdf_output(int num_atoms, arma::irowvec atom_identity, arma::mat coordinates, arma::imat bonding, std::string name, double energy);