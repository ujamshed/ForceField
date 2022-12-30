# ForceField
Implementation of a force field for molecular geometery optimization, based on the MMFF94 by Thomas a. Halgren. Can only be applied to saturated alkanes.

## Installation

### To install ForceField, follow these steps:

Clone the repository using the following command:
```sh
git clone https://github.com/username/ForceField.git
```
Navigate to the directory where you cloned the repository:

```sh
cd ForceField
```

Compile the code using the provided makefile:

```sh
make all
```

## Usage

Use the CL executable ff, found in the bin folder after compilation.

```sh
./ff filename.txt
```

## Creating your own executable

```
#include "../include/input.h"
#include "../include/utility.h"
#include "../include/energy.h"
#include "../include/optimizers.h"

int main()
{
    Molecule Mol(filepath);

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
```

## Dependencies
Forcefield depends on the linear algebra library [Armadillo](https://arma.sourceforge.net/).

## Input Format
View the current input format required by navigating to data/methane_optimal.txt.

1. The first line contains the number of atoms in the molecule
2. The next set of lines are the x, y, z coordinates for each atom and the element symbol.
3. Then comes the bonding information for the atoms (bonding i to j, with bond type) (1=single, 2=double, 3=triple)
4. After bonding information asterisks seperate the structural information including angles, dihedrals, stretch-bend, vanderwaals, etc.

## Updates
1. Calculate structure information interactions including angles, dihedrals, stretch-bending, etc on the fly.
2. Extend the forcefield to all atom types featured in the MMFF94 full implementation.
