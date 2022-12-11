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

    // n1.print("Normal vector 1");
    // n2.print("Normal vector 2");
    // std::cout << "Distance n1: " << dist_n1 << std::endl;
    // std::cout << "Distance n2: " << dist_n2 << std::endl;

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

