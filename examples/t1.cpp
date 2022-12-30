#include "../include/input.h"
#include "../include/utility.h"
#include "../include/energy.h"
#include "../include/optimizers.h"

int main()
{
    std::cout << "Molecule: " << std::endl;
    
    Molecule Ethane("../data/ethane2.txt");
    
    // Get input ethane coordinates
    arma::mat Initial_Ethane_Coordinates = Ethane.getMoleculeCoordinates();

    // Get number of ethane atoms
    int ethane_atoms = Ethane.getNumberOfAtoms();

    // Calculate initial ethane energy
    calcEnergy Ethane_energy(&Ethane);
    std::cout << "Initial Ethane Energy: " << std::endl;
    std::cout << Ethane_energy.total(Initial_Ethane_Coordinates, false) << std::endl;

    // Export initial ethane coordinates (to compare with output)
    Ethane_energy.export_sdf(Initial_Ethane_Coordinates, "ethane_input");

    // Optimize structure and output final energy and export the final output
    //arma::mat output = Ethane_energy.BFGS(0.0001);
    arma::mat output = BFGS(Initial_Ethane_Coordinates, 0.0001, ethane_atoms, Ethane_energy);
    Ethane_energy.total(output);
    Ethane_energy.export_sdf(output, "ethane_out");
}