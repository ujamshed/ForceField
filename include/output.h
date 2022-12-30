#pragma once
#include <armadillo>
#include <fstream>
#include <iomanip>

// Helper function to output the optimized coordinates into an sdf file to view in pymol or related molecule viewer.
void sdf_output(int num_atoms, arma::irowvec atom_identity, arma::mat coordinates, arma::imat bonding, std::string name, double energy);
