#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <armadillo>
#include <fstream>

class Molecule
{
    public:
        int _num_atoms;
        arma::mat _mol_data;
        arma::irowvec _atom_identity;
        arma::irowvec _atom_types;

        // contains bonding information for the molecule (i, j, bond type)
        arma::imat _bonding_data;
        // contains angles information for the molecule (i, j, k)
        arma::imat _angle_data;
        // contains stretch bend information for the molecule (i, j, k)
        arma::imat _stretch_bend_data;
        // contains dihedral angles information for the molecule (i, j, k, l)
        arma::imat _d_angle_data;

        // Maps element symbol to atomic number
        std::map<char, int> symbol_to_int = {{'H', 1}, {'C', 6}};
        
        // Maps atomic number to atom type
        std::map<int, int> an_to_at = {{1, 5}, {6, 1}};

        /**
         * @brief Construct a new Molecule object
         * 
         * @details Holds all the molecule information (atom types, bonds, angles, etc).
         * 
         * @param filename (string): Name of the file which you would like to parse into a Molecule object.
         */
        Molecule(std::string filename);

        /**
         * @brief Get the Molecule Coordinates
         * 
         * @return arma::mat 
         */
        arma::mat getMoleculeCoordinates();

        /**
         * @brief Get the Number Of Atoms in the molecule.
         * 
         * @return int 
         */
        int getNumberOfAtoms();
};
