#include "../include/input.h"
#include "../include/utility.h"
#include "../include/energy.h"

#define IS_TRUE(x) { if ((!x)) std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; }

bool input_coordinates(arma::mat test_coordinates, arma::mat true_coordinates)
{
    return arma::approx_equal(test_coordinates, true_coordinates, "absdiff", 1e-3);
}

bool energy(double test_energy, double true_energy)
{
    return (abs(test_energy-true_energy) < 0.01);
}

int main()
{    
    /*
    * Unit tests on methane coordiantes and initial energy 
    */

    Molecule Methane("../data/methane2.txt");
    
    // Get input ethane coordinates
    arma::mat Initial_methane_Coordinates = Methane.getMoleculeCoordinates();

    arma::mat True_methane_coordinates = {{0.8643,   0.0433,  -0.0827},
                                        {2.0565,   0.0433,  -0.0827},
                                        {0.6002,  -0.7337,  -0.4224},
                                        {0.5002,   0.9491,  -0.5725},
                                        {0.5002,   0.0145,   0.9466}};

    // Check if input coordinates are correct
    IS_TRUE(input_coordinates(Initial_methane_Coordinates, True_methane_coordinates));

    // Calculate initial ethane energy
    calcEnergy Methane_energy(&Methane);
    double meth_energy = Methane_energy.total(Initial_methane_Coordinates, false);

    double true_meth_energy = 23.9428;

    // Check if initial energy is correct
    IS_TRUE(energy(meth_energy, true_meth_energy));


    /*
    * Unit tests on ethane coordiantes and initial energy 
    */

    Molecule ethane("../data/ethane2.txt");
    
    // Get input ethane coordinates
    arma::mat Initial_ethane_Coordinates = ethane.getMoleculeCoordinates();

    arma::mat True_ethane_coordinates = {{1.1605,  -0.1000,  -0.0199},
   {2.7726,  -0.0700,  -0.0199},
   {0.7762,   0.9435,   0.0202},
   {0.8762,  -0.6266,   0.9465},
   {0.7762,  -0.4571,  -0.8264},
   {2.8569,   0.3970,   0.8866},
   {2.9569,   0.2665,  -0.9863},
   {2.9569,  -1.2036,  -0.0600}};

    // Check if input coordinates are correct
    IS_TRUE(input_coordinates(Initial_ethane_Coordinates, True_ethane_coordinates));

    // Calculate initial ethane energy
    calcEnergy ethane_energy(&ethane);
    double eth_energy = ethane_energy.total(Initial_ethane_Coordinates, false);

    double true_eth_energy = 22.401;

    // Check if initial energy is correct
    IS_TRUE(energy(eth_energy, true_eth_energy));
}