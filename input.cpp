#include "input.h"

Molecule::Molecule(std::string filename)
{
    {   
        std::ifstream infile(filename);

        if(!infile.is_open())
        {
            throw std::runtime_error("File does not exist.");
        }

        infile >> _num_atoms;
        std::cout << "Number of atoms is: " << _num_atoms << std::endl;

        _mol_data.zeros(_num_atoms, 3);
        _atom_identity.zeros(_num_atoms);

        std::string line;

        int current_atom = 0;
        char atom_iden;

        // Parsing the atom coordinates in the molecule;
        for (int i = 0; i < _num_atoms; i++)
        {
                infile >> _mol_data(current_atom, 0) >> _mol_data(current_atom, 1) >> _mol_data(current_atom, 2) >> atom_iden;
                
                // Adding atom identity
                //_mol_data(current_atom, 0) = current_atom;

                _atom_identity(current_atom) = symbol_to_int[atom_iden];
                current_atom += 1;
        }

        // Parsing the bonding information in the molecule;
        while (std::getline(infile, line))
        {
            if (line.find("*") != std::string::npos)
            {
                break;
            }

            if (line.size() != 0){
                std::stringstream lineStream(line); 
                
                arma::irowvec bonds;
                bonds.zeros(3);

                lineStream >> bonds(0) >> bonds(1) >> bonds(2);
                _bonding_data = std::move(arma::join_cols(_bonding_data, bonds));
            }
        } 

        // Parsing the angle information in the molecule;
        while (std::getline(infile, line))
        {
            if (line.find("*") != std::string::npos)
            {
                break;
            }

            if (line.size() != 0){
                std::stringstream lineStream(line); 
                
                arma::irowvec angles;
                angles.zeros(3);

                lineStream >> angles(0) >> angles(1) >> angles(2);
                _angle_data = std::move(arma::join_cols(_angle_data, angles));
            }
        }

        // Assigning atom types (trivial implementation assuming only alkanes)
        _atom_types.zeros(_num_atoms);
        for (int i = 0; i < _num_atoms; i++)
        {
            _atom_types(i) = an_to_at[_atom_identity[i]];
        }

        _mol_data.print("Matrix of input data");
        _atom_identity.print("Vector of mol identity");
        _bonding_data.print("Matrix of bonding information");
        _angle_data.print("Matrix of angle information");
        _atom_types.print("Atom Types");
    }
};
