#include "energy.h"

calcEnergy::calcEnergy(Molecule* Mol)
{
    std::cout << "Molecule Loaded!" << std::endl;
    _Mol = Mol;
};

double calcEnergy::bondStretching(double kb, double delta_R)
{
    // cs is the cubic strech constant at -2 A^-1
    double cs = -2;
    return 143.9325 * (kb / 2) * delta_R * delta_R * (1 + cs*delta_R + 7/12*cs*cs*delta_R*delta_R);
}

double calcEnergy::angleBending(double ka, double delta_Angle)
{
    // cb is the cubic bend constant at -0.007 deg^-1 or -0.4 rad^-1
    double cb = -0.007;
    return 0.043844 * (ka/2) * delta_Angle * delta_Angle * (1 + cb* delta_Angle);
}

double calcEnergy::angleStretchBending(double kba_ijk, double kba_kji, double delta_Rij, double delta_Rkj, double delta_Angle)
{
    return 2.51210*(kba_ijk*delta_Rij + kba_kji*delta_Rkj)*delta_Angle;
}

double calcEnergy::outOfPlaneBending(double koop, double X)
{
    return 0.043844 * (koop / 2) * X * X;
}

double calcEnergy::torsion(double V1, double V2, double V3, double omega)
{
    return 0.5* (V1 * (1 + cos(omega)) + V2 * (1 - cos(2*omega)) + V3 * (1 + cos(3*omega)));
}

double calcEnergy::vdw(double dist_Rij, double Ai, double Aj, double alphai, double alphaj, double Ni, double Nj, double Gi, double Gj)
{
    double Rii = Ai * pow(alphai, 0.25);
    double Rjj = Aj * pow(alphaj, 0.25);
    double yij = (Rii - Rjj) / (Rii + Rjj);
    
    double Rij = 0.5 * (Rii + Rjj) * (1 + 0.2*(1-exp(-12*pow(yij, 2))));
    double epij = ((181.16 * Gi * Gj * alphai * alphaj) / (pow((alphai / Ni), 0.5) + pow((alphaj / Nj), 0.5))) * (1 / pow(Rij, 6));

    return epij * pow((1.07*Rij)/(dist_Rij + 0.07*Rij) , 7) * ((1.12*pow(Rij, 7) / (pow(dist_Rij, 7) + 0.12*pow(Rij, 7))) - 2);
}

double calcEnergy::electrostatic(double delta_R, double qi, double qj)
{
    // dielectric constant of 1 is chosen to assume its within vacuum. From introduction to computational chemistry page 41.
    double dielectric_constant = 1.0;
    double delta = 0.05;

    return 332.0716*qi*qj / (dielectric_constant * (delta_R + delta));
}

double calcEnergy::bondingContributions(arma::mat molecule_coordinates, bool verbose)
{
    // Loop through all the bonding data and evaluate the bonding energy summation
    double total = 0;

    for (int i=0; i < _Mol->_bonding_data.n_rows; i++)
    {
        // Get indices for the atoms involved in the bonding
        int atom_i_index = _Mol->_bonding_data(i, 0);
        int atom_j_index = _Mol->_bonding_data(i, 1);

        // Get coordinates for the atoms involved in the bonding
        arma::rowvec atom_i_coordinates = molecule_coordinates.row(atom_i_index);
        arma::rowvec atom_j_coordinates = molecule_coordinates.row(atom_j_index);

        // Get atom types for the atoms involved in the bonding
        int atom_i_type = _Mol->_atom_types(atom_i_index);
        int atom_j_type = _Mol->_atom_types(atom_j_index);

        // Get distance between the 2 atoms
        double dist = euclidian_distance(atom_i_coordinates, atom_j_coordinates);
    
        // Get params
        std::pair<int, int> atom_types = {atom_i_type, atom_j_type};
        std::pair<double, double> params = bs_params[atom_types];

        if (verbose)
        {
            std::cout << "Atoms: " << atom_i_index << " " << atom_j_index << " " << " Distance: " << dist << std::endl;
        }

        double delta_R = dist - params.second;
        total += bondStretching(params.first, delta_R);
    }
    //std::cout << "Bonding Total: " << total << std::endl;
    return total;
}

double calcEnergy::angleContributions(arma::mat molecule_coordinates, bool verbose)
{
    // Loop through all the angle data and evaluate the bonding energy summation
    double total = 0;

    for (int i=0; i < _Mol->_angle_data.n_rows; i++)
    {
        // Get indices for the atoms involved in the angle
        int atom_i_index = _Mol->_angle_data(i, 0);
        int atom_j_index = _Mol->_angle_data(i, 1);
        int atom_k_index = _Mol->_angle_data(i, 2);

        // Get coordinates for the atoms involved in the bonding
        arma::rowvec atom_i_coordinates = molecule_coordinates.row(atom_i_index);
        arma::rowvec atom_j_coordinates = molecule_coordinates.row(atom_j_index);
        arma::rowvec atom_k_coordinates = molecule_coordinates.row(atom_k_index);

        // Get atom types for the atoms involved in the bonding
        int atom_i_type = _Mol->_atom_types(atom_i_index);
        int atom_j_type = _Mol->_atom_types(atom_j_index);
        int atom_k_type = _Mol->_atom_types(atom_k_index);

        // Get angle between the 3 atoms
        double angle = calc_angle(atom_i_coordinates, atom_j_coordinates, atom_k_coordinates);
    
        // Get params
        std::tuple<int, int, int> atom_types = {atom_i_type, atom_j_type, atom_k_type};
        std::pair<double, double> params = ab_params[atom_types];


        if (verbose)
        {
            std::cout << "Atoms: " << atom_i_index << " " << atom_j_index << " " << atom_k_index << " Angle: " << angle << std::endl;
        }

        double delta_Angle = angle - params.second;
        total += angleBending(params.first, delta_Angle);
    }
    //std::cout << "Angle Bending Total: " << total << std::endl;
    return total;
}

double calcEnergy::angleStretchBendingContributions(arma::mat molecule_coordinates)
{
    // Loop through all the angle data and evaluate the bonding energy summation
    double total = 0;

    for (int i=0; i < _Mol->_stretch_bend_data.n_rows; i++)
    {
        // Get indices for the atoms involved in the angle
        int atom_i_index = _Mol->_stretch_bend_data(i, 0);
        int atom_j_index = _Mol->_stretch_bend_data(i, 1);
        int atom_k_index = _Mol->_stretch_bend_data(i, 2);

        // Get coordinates for the atoms involved in the bonding
        arma::rowvec atom_i_coordinates = molecule_coordinates.row(atom_i_index);
        arma::rowvec atom_j_coordinates = molecule_coordinates.row(atom_j_index);
        arma::rowvec atom_k_coordinates = molecule_coordinates.row(atom_k_index);

        // Get atom types for the atoms involved in the bonding
        int atom_i_type = _Mol->_atom_types(atom_i_index);
        int atom_j_type = _Mol->_atom_types(atom_j_index);
        int atom_k_type = _Mol->_atom_types(atom_k_index);

        // Get distance between atoms i & j and atoms j & k
        double R_ij = euclidian_distance(atom_i_coordinates, atom_j_coordinates);
        double R_jk = euclidian_distance(atom_j_coordinates, atom_k_coordinates);

        // Get angle between the 3 atoms
        double angle = calc_angle(atom_i_coordinates, atom_j_coordinates, atom_k_coordinates);

        // Get bond stretching params
        std::pair<int, int> ij_atom_types = {atom_i_type, atom_j_type};
        std::pair<double, double> bs_ij_params = bs_params[ij_atom_types];

        std::pair<int, int> jk_atom_types = {atom_j_type, atom_k_type};
        std::pair<double, double> bs_jk_params = bs_params[jk_atom_types];

        // Get angle params
        std::tuple<int, int, int> ijk_atom_types = {atom_i_type, atom_j_type, atom_k_type};
        std::pair<double, double> ab = ab_params[ijk_atom_types];
    
        // Get angle strech bend params
        std::pair<double, double> asb = asb_params[ijk_atom_types];

        double delta_Angle = (angle - ab.second);
        double delta_Rij = R_ij - bs_ij_params.second;
        double delta_Rjk = R_jk - bs_jk_params.second;

        total += angleStretchBending(asb.first, asb.second, delta_Rij, delta_Rjk, delta_Angle);
    }

    //std::cout << "Angle Stretch Bending Total: " << total << std::endl;
    return total;
}

// Saturated aliphatics do not have oop contributions, so this is not implemented for now.
void calcEnergy::oopContributions()
{

};

double calcEnergy::torsionalContributions(arma::mat molecule_coordinates)
{
    // Loop through all the angle data and evaluate the bonding energy summation
    double total = 0;

    for (int i=0; i < _Mol->_d_angle_data.n_rows; i++)
    {
        // Get indices for the atoms involved in the torsional angle
        int atom_i_index = _Mol->_d_angle_data(i, 0);
        int atom_j_index = _Mol->_d_angle_data(i, 1);
        int atom_k_index = _Mol->_d_angle_data(i, 2);
        int atom_l_index = _Mol->_d_angle_data(i, 3);

        // Get coordinates for the atoms involved in the bonding
        arma::rowvec atom_i_coordinates = molecule_coordinates.row(atom_i_index);
        arma::rowvec atom_j_coordinates = molecule_coordinates.row(atom_j_index);
        arma::rowvec atom_k_coordinates = molecule_coordinates.row(atom_k_index);
        arma::rowvec atom_l_coordinates = molecule_coordinates.row(atom_l_index);

        // Get atom types for the atoms involved in the bonding
        int atom_i_type = _Mol->_atom_types(atom_i_index);
        int atom_j_type = _Mol->_atom_types(atom_j_index);
        int atom_k_type = _Mol->_atom_types(atom_k_index);
        int atom_l_type = _Mol->_atom_types(atom_l_index);

        // Get angle between the 4 atoms
        double angle = calc_dihedral(atom_i_coordinates, atom_j_coordinates, atom_k_coordinates, atom_l_coordinates);

        // Get torsional angle params
        std::tuple<int, int, int, int> ijkl_atom_types = {atom_i_type, atom_j_type, atom_k_type, atom_l_type};
        std::tuple<double, double, double> tp = t_params[ijkl_atom_types];

        // Debugging
        //std::cout << "Atoms: " << atom_i_index << " " << atom_j_index << " " << atom_k_index << " " << atom_l_index << " Angle: " << angle << std::endl;

        total += torsion(std::get<0>(tp), std::get<1>(tp), std::get<2>(tp), angle);
    }

    //std::cout << "Torsion Total: " << total << std::endl;
    return total;
}

double calcEnergy::vdwContributions(arma::mat molecule_coordinates)
{
    // Loop through all the dihedral angles (because they are at a minimum 3 bonds away) and get the interactions between the furthest 2
    // atoms (Will only work for ethane).
    double total = 0;

    for (int i=0; i < _Mol->_d_angle_data.n_rows; i++)
    {

        // Get index for atom i involved in the interaction (the first atom in the bonding)
        int atom_i_index = _Mol->_d_angle_data(i, 0);

        // Get coordinates for atom i involved in the interaction
        arma::rowvec atom_i_coordinates = molecule_coordinates.row(atom_i_index);

        // Get atom type for atom i involved in the interaction
        int atom_i_type = _Mol->_atom_types(atom_i_index);

        // Get indices for atom j involved in the interaction
        int atom_j_index = _Mol->_d_angle_data(i, 3);

        // Get coordinates for atom j involved in the interaction
        arma::rowvec atom_j_coordinates = molecule_coordinates.row(atom_j_index);

        // Get atom type for atom j involved in the interaction
        int atom_j_type = _Mol->_atom_types(atom_j_index);

        // Get distance between the 2 atoms
        double dist_Rij = euclidian_distance(atom_i_coordinates, atom_j_coordinates);

        // Debugging
        // atom_i_coordinates.print("Atom i");
        // std::cout << "Atom i index: " << atom_i_index << " atom i type: " << atom_i_type << std::endl;
        // atom_j_coordinates.print("Atom j");
        // std::cout << "Atom j index: " << atom_j_index << " atom j type: " << atom_j_type << std::endl;
        // std::cout << std::endl;

        // Get params
        std::tuple<double, double, double, double> i_params = vdw_params[atom_i_type];
        std::tuple<double, double, double, double> j_params = vdw_params[atom_j_type];

        double alpha_i = std::get<0>(i_params);
        double alpha_j = std::get<0>(j_params);
        double N_i = std::get<1>(i_params);
        double N_j = std::get<1>(j_params);
        double Ai = std::get<2>(i_params);
        double Aj = std::get<2>(j_params);
        double Gi = std::get<3>(i_params);
        double Gj = std::get<3>(j_params);

        total += vdw(dist_Rij, Ai, Aj, alpha_i, alpha_j, N_i, N_j, Gi, Gj);             
    }
    //std::cout << "Vdw Total: " << total << std::endl;
    return total;
};

// Saturated aliphatics do not have electrostatic contributions, so this is not implemented for now.
void calcEnergy::electrostaticContributions()
{
};

double calcEnergy::total(arma::mat molecule_coordinates, bool verbose)
{
    double total = bondingContributions(molecule_coordinates, verbose) + angleContributions(molecule_coordinates, verbose) + angleStretchBendingContributions(molecule_coordinates) + torsionalContributions(molecule_coordinates) + vdwContributions(molecule_coordinates);
    //std::cout << "Total is: " << total << std::endl;
    return total;
}

void calcEnergy::export_sdf(arma::mat coordinates, std::string name)
{
    double energy = total(coordinates);
    sdf_output(_Mol->_num_atoms, _Mol->_atom_identity, coordinates, _Mol->_bonding_data, name, energy);
}
