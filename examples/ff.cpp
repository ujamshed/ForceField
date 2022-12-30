#include "../include/input.h"
#include "../include/utility.h"
#include "../include/energy.h"
#include "../include/optimizers.h"
#include <stdexcept>

int main(int argc, const char* argv[])
{
    if (argc < 2)
    {
        throw std::invalid_argument( "Please enter the file path for the molecule you wish to optimize" );
    }

    std::cout << "Molecule: " << std::endl;
    
    Molecule Mol(argv[1]);
    
    // Get input coordinates
    arma::mat Initial_Mol_Coordinates = Mol.getMoleculeCoordinates();

    // Get atoms
    int atoms = Mol.getNumberOfAtoms();

    // Calculate initial energy
    calcEnergy Mol_energy(&Mol);
    std::cout << "Initial Mol Energy: " << std::endl;
    std::cout << Mol_energy.total(Initial_Mol_Coordinates, false) << std::endl;

    // Export initial coordinates (to compare with output)
    Mol_energy.export_sdf(Initial_Mol_Coordinates, "mol_input");

    // Optimize structure and output final energy and export the final output
    //arma::mat output = Mol_energy.BFGS(0.0001);
    arma::mat output = BFGS(Initial_Mol_Coordinates, 0.0001, atoms, Mol_energy);
    
    Mol_energy.total(output);
    Mol_energy.export_sdf(output, "mol_out");
}