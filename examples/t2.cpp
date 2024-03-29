#include "../include/input.h"
#include "../include/utility.h"
#include "../include/energy.h"
#include "../include/optimizers.h"

int main()
{
    std::cout << "Molecule: " << std::endl;
    
    Molecule methane("../data/methane2.txt");
    
    // Get input methane coordinates
    arma::mat Initial_Methane_Coordinates = methane.getMoleculeCoordinates();

    // Get methane atoms
    int methane_atoms = methane.getNumberOfAtoms();

    // Calculate initial methane energy
    calcEnergy Methane_energy(&methane);
    std::cout << "Initial Methane Energy: " << std::endl;
    std::cout << Methane_energy.total(Initial_Methane_Coordinates, false) << std::endl;

    // Export initial methane coordinates (to compare with output)
    Methane_energy.export_sdf(Initial_Methane_Coordinates, "methane2_input");

    // Optimize structure and output final energy and export the final output
    arma::mat output = BFGS(Initial_Methane_Coordinates, 0.0001, methane_atoms, Methane_energy);
    
    Methane_energy.total(output);
    Methane_energy.export_sdf(output, "Methane2_out");
}