#include "utility.h"

arma::rowvec vector_subtraction(arma::rowvec atom_i, arma::rowvec atom_j)
{
    arma::rowvec difference = {atom_i[0] - atom_j[0], atom_i[1] - atom_j[1], atom_i[2] - atom_j[2]};
    return difference;
};

double euclidian_distance(arma::rowvec atom_i, arma::rowvec atom_j)
{
    return sqrt(pow((atom_i[0] - atom_j[0]), 2) + pow((atom_i[1] - atom_j[1]), 2) + pow((atom_i[2] - atom_j[2]), 2));
};

double calc_angle(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k)
{
    double d12 = euclidian_distance(atom_i, atom_j);
    double d23 = euclidian_distance(atom_j, atom_k);
    double d13 = euclidian_distance(atom_i, atom_k);

    double angle = acos((-pow(d13, 2) + pow(d12, 2) + pow(d23, 2)) / (2 * d12 * d23));

    return angle * (180/M_PI);
};

double calc_dihedral(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k, arma::rowvec atom_l)
{

    // Get vectors joining the 4 points together
    arma::rowvec v1 = vector_subtraction(atom_i, atom_j);
    arma::rowvec v2 = vector_subtraction(atom_j, atom_k);
    arma::rowvec v3 = vector_subtraction(atom_k, atom_l);

    // Get normal vectors
    arma::rowvec n1 = arma::cross(v1, v2);
    arma::rowvec n2 = arma::cross(v2, v3);

    // Distance from origin to normal vectors to get their length
    arma::rowvec origin = {0, 0, 0};

    double dist_n1 = euclidian_distance(n1, origin);
    double dist_n2 = euclidian_distance(n2, origin);

    double value = arma::dot(n1, n2) / (dist_n1 * dist_n2);

    // Check to see if the value is less than -1 or greater than 1 (limits for acos function)
    // if (value < -1.00)
    // {
    //     value = -1.00;
    // }
    // if (value > 1.00)
    // {
    //     value = 1.00;
    // }

    double angle = acos(value);

    return angle * (180/M_PI);
};

double calc_partial_atomic_charge(arma::mat atom_types)
{
    // Initializing q0 at 0.0 as per the paper.
    double q0 = 0.0;
    return q0;
};

void sdf_output(int num_atoms, arma::irowvec atom_identity, arma::mat coordinates, arma::imat bonding)
{
    std::ofstream ofile("output.sdf");

    std::map<int, char> int_to_ele_char = {{1, 'H'}, {6, 'C'}};

    if (ofile.is_open())
    {
        // Add necessary spaces and number of atoms and number of bonds
        ofile << "\n";
        ofile << "\n";
        ofile << "\n";
        ofile << " " << num_atoms << " " << bonding.n_rows << "\n";
        
        // Get all the coordinate information
        for (int i=0; i < num_atoms; i++)
        {
            ofile << " " << coordinates(i, 0) << " " << coordinates(i, 1) << " " << coordinates(i, 2) << " " << int_to_ele_char[atom_identity(i)] << std::endl;
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