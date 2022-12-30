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
    double dij = euclidian_distance(atom_i, atom_j);
    double djk = euclidian_distance(atom_j, atom_k);
    double dik = euclidian_distance(atom_i, atom_k);

    double angle = acos((-pow(dik, 2) + pow(dij, 2) + pow(djk, 2)) / (2 * dij * djk));

    return angle * (180/M_PI);
};

double calc_dihedral(arma::rowvec atom_i, arma::rowvec atom_j, arma::rowvec atom_k, arma::rowvec atom_l)
{

    // Calculating dihedral angle as per eq 12.8.5 in Crystal Structure Analysis for Chemists and Biologists

    // Get distances between all the coordinates
    double d12 = euclidian_distance(atom_i, atom_j);
    double d13 = euclidian_distance(atom_i, atom_k); 
    double d14 = euclidian_distance(atom_i, atom_l); 
    double d23 = euclidian_distance(atom_j, atom_k);
    double d24 = euclidian_distance(atom_j, atom_l);
    double d34 = euclidian_distance(atom_k, atom_l);

    // Calculate P
    double P = pow(d12, 2) * (pow(d23, 2) + pow(d34,2) - pow(d24, 2)) + pow(d23,2)*(-pow(d23,2) + pow(d34,2) + pow(d24,2)) + pow(d13, 2)*(pow(d23,2) - pow(d34,2) + pow(d24, 2)) - 2*pow(d23, 2)*pow(d14,2);

    // Calculate Q
    double Q = ((d12 + d23 + d13) * (d12 + d23 - d13) * (d12 - d23 + d13) * (-d12 + d23 + d13) * (d23 + d34 + d24) * (d23 + d34 -d24) * (d23 - d34 + d24) * (-d23 + d34 + d24));

    double value = P/sqrt(Q);
    if ((P/sqrt(Q)) < -1)
    {
        value = -1;
    }
    if ((P/sqrt(Q)) > 1)
    {
        value = 1;
    }
    return acos(value);
};

double calc_partial_atomic_charge(arma::mat atom_types)
{
    // Initializing q0 at 0.0 as per the paper.
    double q0 = 0.0;
    return q0;
};
