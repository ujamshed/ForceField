#include "output.h"

void sdf_output(int num_atoms, arma::irowvec atom_identity, arma::mat coordinates, arma::imat bonding, std::string name, double energy)
{
    std::string filename = "../output/" + name + ".sdf";
    std::ofstream ofile(filename);

    std::map<int, char> int_to_ele_char = {{1, 'H'}, {6, 'C'}};

    if (ofile.is_open())
    {
        // Add necessary spaces and number of atoms and number of bonds
        ofile << "MMFF94 Implementation as CHEM279 final project\n";
        ofile << "Configuration Energy: " << energy << "\n";
        ofile << "\n";
        ofile << " " << num_atoms << " " << bonding.n_rows << "\n";
                
        // Get all the coordinate information
        for (int i=0; i < num_atoms; i++)
        {
            for (int j=0; j<3; j++)
            {
                if (coordinates(i,j) < 0)
                {
                    // 3 spaces if the number is negative
                    ofile << std::fixed << std::setprecision(4) << "   " << coordinates(i, j);
                }
                else
                {
                    // 4 spaces otherwise
                    ofile << std::fixed << std::setprecision(4) << "    " << coordinates(i, j);
                }
            }
            ofile << " " << int_to_ele_char[atom_identity(i)] << std::endl;
        }

        // Get the bonding information
        for (int j=0; j < bonding.n_rows; j++)
        {
            // Need to increment by 1 for the atoms since atoms are numbered starting at 1 in sdf file format.
            ofile << "  " << bonding(j, 0)+1 << "  " << bonding(j, 1)+1 << "  " << bonding(j, 2) << std::endl;
        }
        ofile.close();
    }
}