#include "utility.h"

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