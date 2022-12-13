#include "input.h"
#include "utility.h"

class calcEnergy
{
    private:
        Molecule* _Mol;

        // Maps atom types to bond stretching parameters
        // First pair is the atom types involved in the bond, and the second pair are the params
        // First element of the second pair is the kb, while the second element of the second pair is the r0
        // Need 1,5 and 5,1 because C-H and H-C when considering the stretch bending contributions.
        std::map<std::pair<int, int>, std::pair<double, double>> bs_params = { {{1,1} , {4.258, 1.508}}, {{1,5} , {4.766, 1.093}}, {{5,1} , {4.766, 1.093}} };
        
        // Maps atom types to angle parameters
        // First tuple is the atom types involved in the angle, and the second pair are the params
        // First element of the second pair is the ka, while the second element of the second pair is the a0
        std::map<std::tuple<int, int, int>, std::pair<double, double>> ab_params = { {{1,1,5} , {0.636, 110.549}}, {{5, 1, 5}, {0.516, 108.836}} };
        
        // Maps atom types to stretch bond parameters
        // First tuple is the atom types involved in the interaction and second pair are the params
        // First element of the second pair is the kbaIJK, while the second element of the second pair is the kbaKJI
        std::map<std::tuple<int, int, int>, std::pair<double, double>> asb_params = { {{1,1,5} , {0.227, 0.070}}, {{5, 1, 5}, {0.115, 0.115}} };

        // Maps atom types to torsion parameters
        // First tuple is the atom types involved in the interaction and the second pair are the params
        // Elements of the second tuple are V1, V2, V3 respectively
        std::map<std::tuple<int, int, int, int>, std::tuple<double, double, double>> t_params = { {{5,1,1,5} , {0.284, -1.386, 0.314}} };

        // Maps atom types to vdw parameters
        // Int is the atom type involved in the interaction
        // Elements of the second tuple are alpha, N, A and G, respectively
        std::map<int, std::tuple<double, double, double, double>> vdw_params = { {1, {1.050, 2.490, 3.890, 1.282}}, {5, {0.250, 0.800, 4.200, 1.209}} };


    public:
        
        calcEnergy(Molecule* Mol)
        {
            std::cout << "Molecule Loaded!" << std::endl;
            _Mol = Mol;
        };

        // Need helper function to calculate the distance between the two atoms i and j.
        // Calculates the bond stretching energy between 2 atoms, i and j.
        double bondStretching(double kb, double delta_R)
        {
            // cs is the cubic strech constant at -2 A^-1
            double cs = -2;
            return 143.9325 * (kb / 2) * delta_R * delta_R * (1 + cs*delta_R + 7/12*cs*cs*delta_R*delta_R);
        }

        // Need helper function to calculate angle between the 3 atoms
        // delta_Angle is angle_ijk - reference angle.
        // Calculates the angle bending between atoms i, j and k, where j is the middle atom.
        double angleBending(double ka, double delta_Angle)
        {
            // cb is the cubic bend constant at -0.007 deg^-1 or -0.4 rad^-1
            double cb = -0.007;
            return 0.043844 * (ka/2) * delta_Angle * delta_Angle * (1 + cb* delta_Angle);
        }

        // Calculates the angle stretch bending between atoms i, j, k. Distances between i, j and k, j need to be calculated, as well
        // as the angle between i, j, k.
        double angleStretchBending(double kba_ijk, double kba_kji, double delta_Rij, double delta_Rkj, double delta_Angle)
        {
            return 2.51210*(kba_ijk*delta_Rij + kba_kji*delta_Rkj)*delta_Angle;
        }

        // Need helper function to calculate the wilson angle (X) between the bond j-l, and the plane i-j-k.
        // Calculates the out of plane bending at atom j in i,j,k.
        double outOfPlaneBending(double koop, double X)
        {
            return 0.043844 * (koop / 2) * X * X;
        }

        // Need helper function to calculate the torsional angle (omega) between atoms i, j, k, l.
        // Calculates the torsional energy based on atoms i, j, k, l.
        double torsion(double V1, double V2, double V3, double omega)
        {
            return 0.5* (V1 * (1 + cos(omega)) + V2 * (1 - cos(2*omega)) + V3 * (1 + cos(3*omega)));
        }

        // Calculates the vanderwaals interaction based on the buffered-14-7 form, based on atoms i and j.
        double vdw(double dist_Rij, double Ai, double Aj, double alphai, double alphaj, double Ni, double Nj, double Gi, double Gj)
        {
            double Rii = Ai * pow(alphai, 0.25);
            double Rjj = Aj * pow(alphaj, 0.25);
            double yij = (Rii - Rjj) / (Rii + Rjj);
            
            double Rij = 0.5 * (Rii + Rjj) * (1 + 0.2*(1-exp(-12*pow(yij, 2))));
            double epij = ((181.16 * Gi * Gj * alphai * alphaj) / (pow((alphai / Ni), 0.5) + pow((alphaj / Nj), 0.5))) * (1 / pow(Rij, 6));

            return epij * pow((1.07*Rij)/(dist_Rij + 0.07*Rij) , 7) * ((1.12*pow(Rij, 7) / (pow(dist_Rij, 7) + 0.12*pow(Rij, 7))) - 2);
        }

        // Need helper function to calculate the partial atomic charges (qi, qj) based on bond charge increments (bci)
        // Calculates the electrostatic interaction between atoms i, j.
        double electrostatic(double delta_R, double qi, double qj)
        {
            // dielectric constant of 1 is chosen to assume a vacuum. From introduction to computational chemistry page 41.
            double dielectric_constant = 1.0;
            double delta = 0.05;

            return 332.0716*qi*qj / (dielectric_constant * (delta_R + delta));
        }

        double bondingContributions(arma::mat molecule_coordinates, bool verbose=false)
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

        double angleContributions(arma::mat molecule_coordinates, bool verbose = false)
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

        double angleStretchBendingContributions(arma::mat molecule_coordinates)
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
        
        // Ethane does not have oop contributions, so this is not implemented for now.
        void oopContributions()
        {

        };

        double torsionalContributions(arma::mat molecule_coordinates)
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

                //std::cout << "Atom i: " << std::get<0>(ijkl_atom_types) << " Atom j: " << std::get<1>(ijkl_atom_types) << " Atom k: " << std::get<2>(ijkl_atom_types) << " Atom l: " << std::get<3>(ijkl_atom_types) << std::endl;

                total += torsion(std::get<0>(tp), std::get<1>(tp), std::get<2>(tp), angle);
            }

            //std::cout << "Torsion Total: " << total << std::endl;
            return total;
        }

        double vdwContributions(arma::mat molecule_coordinates)
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

        // Ethane does not have electrostatic contributions, so this is not implemented for now.
        void electrostaticContributions()
        {
        };

        double total(arma::mat molecule_coordinates, bool verbose=false)
        {
            double total = bondingContributions(molecule_coordinates, verbose) + angleContributions(molecule_coordinates, verbose) + angleStretchBendingContributions(molecule_coordinates) + torsionalContributions(molecule_coordinates) + vdwContributions(molecule_coordinates);
            std::cout << "Total is: " << total << std::endl;
            return total;
        }

        // Function to calculate central difference using input matrix
        arma::mat finite_central_differences_sd(arma::mat Coordinates, double step_size)
        {
            arma::mat zero;
            zero.zeros(_Mol->_num_atoms,3);

            // Loop through all atoms.
            for (int i=0; i < _Mol->_num_atoms; i++)
            {
                // Loop through x y and z coordinates to get the gradient at each one keeping the rest constant.
                for (int k=0; k < 3; k++)
                {
                    // Get current atom coordinates (x, y, z).
                    arma::mat new_coordinates_plus = Coordinates;
                    arma::mat new_coordinates_minus = Coordinates;

                    // Step size adjustments for x, y or z.
                    new_coordinates_plus(i, k) += step_size;
                    new_coordinates_minus(i, k) -= step_size;

                    // Calculate the plus and minus energies based on the movement of the atoms
                    double plus_E = total(new_coordinates_plus);
                    double minus_E = total(new_coordinates_minus);

                    double force = -(plus_E - minus_E) / (2*step_size);
                    zero(i, k) = force;
                }
            }
            zero.print("Current Gradient");
            return zero;
        }

        // Steepest Descent
        arma::mat steepest_descent(double step_size, double tol)
        {
            arma::mat initial_guess = _Mol->_mol_data;

            initial_guess.print("Initial Guess");

            arma::mat derivative = finite_central_differences_sd(initial_guess, step_size);
            double initial_energy = total(initial_guess);
            std::cout << "Initial Energy: " << initial_energy << std::endl;

            int count = 0;

            
            while (arma::norm(derivative, 2) > tol)
            {
                arma::mat new_positions = initial_guess + (derivative / arma::norm(derivative, 2)) * step_size;
            
                if (total(new_positions) < total(initial_guess))
                {
                    initial_guess = new_positions;
                    step_size *= 1.2;
                    derivative = finite_central_differences_sd(initial_guess, step_size);

                    std::cout << std::endl;
                    std::cout << "Step: " << count << std::endl;
                    std::cout << "Energy: " << total(initial_guess) << std::endl;
                    std::cout << "Force Norm: " << arma::norm(derivative, 2) << std::endl;
                    initial_guess.print("New positions");
                    std::cout << "New step size: " << step_size << std::endl;
                    std::cout << std::endl;

                }

                else
                {
                    step_size *= 0.5;
                }

            count += 1;
            }
        return initial_guess;
        }

        // Steepest Descent with Line search
        // Needs work
        arma::mat steepest_descent_line_search(double step_size, double tol)
        {
            double gr = (sqrt(5) + 1) / 2;
            arma::mat initial_guess = _Mol->_mol_data;

            initial_guess.print("Initial Guess");

            arma::mat derivative = finite_central_differences_sd(initial_guess, 0.0001);
            double initial_energy = total(initial_guess);
            std::cout << "Initial Energy: " << initial_energy << std::endl;

            int count = 0;

            while (arma::norm(derivative, 2) > tol)
            {
                // Reset step_size so it can be minimized towards tolerance again
                step_size = 1.0;

                // Get new position
                arma::mat new_positions = initial_guess + (derivative / arma::norm(derivative, 2)) * step_size;

                // Bracketing step
                // Bracketing stops when lj potential of new_positions is higher than initial_positions
                while (total(new_positions) < total(initial_guess))
                {
                    step_size *= 1.2;
                    new_positions = initial_guess + (derivative / arma::norm(derivative, 2)) * step_size;
                    std::cout << "Currently Bracketing..." << std::endl;
                }

                // Golden section
                // A = initial guess
                // B = new position
                arma::mat C = new_positions - (new_positions - initial_guess) / gr; 
                arma::mat D = initial_guess + (new_positions - initial_guess) / gr;
                
                // line search to optimize step_size
                while (step_size > tol)
                {
                    if (total(C) < total(D))
                    {
                        new_positions = D;
                    }
                    
                    else
                    {
                        initial_guess = C;
                    }

                    C = new_positions - (new_positions - initial_guess) / gr; 
                    D = initial_guess + (new_positions - initial_guess) / gr;

                    // solve for step_size, not necessary
                    arma::mat difference = (new_positions - initial_guess);
                    arma::mat unit_vector = (derivative / arma::norm(derivative, 2));

                    step_size = difference(0, 2) / unit_vector(0, 2);
                }
                
                // Update initial guess to the new optimal
                initial_guess = (new_positions + initial_guess) / 2;

                derivative = finite_central_differences_sd(initial_guess, 0.0001);
                std::cout << std::endl;
                std::cout << "Step: " << count << std::endl;
                std::cout << "Energy: " << total(initial_guess) << std::endl;
                std::cout << "Force Norm: " << arma::norm(derivative, 2) << std::endl;
                initial_guess.print("New positions");
                std::cout << std::endl;

                count += 1;
            }
            return initial_guess;
        }

        // Function to export the molecule information in sdf format
        void export_sdf(arma::mat coordinates)
        {
            sdf_output(_Mol->_num_atoms, _Mol->_atom_identity, coordinates, _Mol->_bonding_data);
        }
};

int main()
{
    std::cout << "First Molecule: " << std::endl;
    
    //Molecule Ethane("data/methane.txt");
    Molecule Ethane("data/methane2.txt");
    arma::mat Ethane_Coordinates = Ethane.getMoleculeCoordinates();
    calcEnergy Ethane_energy(&Ethane);
    std::cout << "Initial Methane Energy: " << std::endl;
    Ethane_energy.total(Ethane_Coordinates, true);
    //Ethane_energy.finite_central_differences_sd(Ethane_Coordinates, 0.0001);
    arma::mat output = Ethane_energy.steepest_descent(0.001, 0.01);
    
    //arma::mat output = Ethane_energy.steepest_descent_line_search(0.001, 0.01);
    output.print("Optimized coordinates");
    
    Ethane_energy.total(output, true);
    Ethane_energy.export_sdf(output);
}